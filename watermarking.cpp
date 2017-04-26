// #include "watermarking.h"
#include <vector>
#include <Magick++.h> 
#include <iostream>
#include <cstdlib> //atoi
#include <cmath>

#define BLK_SIZE 8

using namespace std; 
using namespace Magick;

// matrix_t<T> type
template< typename T > using matrix_t = vector< vector<T> >;

// matrix3d_t<T> type
template< typename T > using matrix3d_t = vector< matrix_t<T> >;

// function declarations
template< typename T >
matrix_t<T> get_block(const matrix_t<T> &matrix, int i, int j, int n=BLK_SIZE);
template< typename T >
matrix3d_t<T> get_blocks(const matrix_t<T> &matrix, int n=BLK_SIZE);
template< typename T >
void concat_blocks(matrix_t<T> &result, matrix3d_t<T> &blocks);
Image load(char *filename);
matrix_t<unsigned char> image2matrix(Image &image);
template< typename T>
void matrix2image(const matrix_t<T> &matrix, Image &image);
template< typename T >
void embed(matrix_t<T> &matrix, vector<T> &watermark, int strength, int iteration, int size=BLK_SIZE);
template< typename T >
vector<T> extract(matrix_t<T> &matrix, int strength, int size=BLK_SIZE);
template< typename T >
void mat_init(matrix_t<T> &matrix, int size);
template< typename T >
matrix_t<T> mat_multiply(matrix_t<T> &a, matrix_t<T> &b, int size=BLK_SIZE);
template< typename T >
matrix_t<T> mat_transpose(matrix_t<T> &matrix);
template< typename T >
double dot_product(const vector<T> &vec1, const vector<T> &vec2);
template< typename T >
double normalize(const vector<T> &vec);
template< typename T >
void vector_concat(const vector<T> &vec, vector<T> &result);
template< typename T >
void gram_schmidt(matrix_t<T> &matrix, matrix3d_t<T> &q_list, matrix3d_t<T> &r_list);
template< typename T >
double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm);
template< typename T >
double nc(const vector<T> &orig, const vector<T> &extr);
template< typename T > 
void print(const matrix_t<T> &matrix);


// get n x n block from matrix on i,j position
// throws out_of_range on error
template< typename T >
matrix_t<T> get_block(const matrix_t<T> &matrix, int i, int j, int n/*=BLK_SIZE*/)
{
    matrix_t<T> block(n);

    for(int k = 0; k < n; ++k)
        for(int l = 0; l < n; ++l)
            block[k].emplace_back(matrix.at(k+i).at(l+j));

    return block;
}

// ???
template< typename T >
matrix3d_t<T> get_blocks(const matrix_t<T> &matrix, int n/*=BLK_SIZE*/)
{
    matrix_t<T> block(n);
    matrix3d_t<T> blocks;
    for (int i = 0; i < matrix.size(); i += n)
        for (int j = 0; j < matrix.size(); j += n) {
            block.clear();
            for(int k = 0; k < n; ++k)
                for(int l = 0; l < n; ++l)
                    block[k].emplace_back(matrix.at(k+i).at(l+j));
            blocks.emplace_back(block);
        }

    return blocks;
}

template< typename T >
void concat_blocks(matrix_t<T> &result, matrix3d_t<T> &blocks)
{
    int i = 0;
    int j = 0;
    for(const auto block : blocks) {
        if (j < result.size()) {
            j += BLK_SIZE;
        }
        else {
            j = 0;
            i++;
        }
        for (int k = 0; k < block.size(); ++k)
            for (int l = 0; l < block.size(); ++l)
                result.at(k+i).at(l+j) = block.at(k).at(l);    
    }

    // for (int k = 0; k < block.size(); ++k)
    //     for (int l = 0; l < block.size(); ++l)
    //         result.at(k+i).at(l+j) = block.at(k).at(l);
}

// load input image
Image load(char *filename)
{
    Image image(filename);
    return image;
}

// transform Image objet to matrix_t
// every value is in range <0-255>
//   0     - black
//   255   - white
// template< typename T>
matrix_t<unsigned char> image2matrix(Image &image)
{
    int w = image.columns();
    int h = image.rows();
    int n = w; // image is n x n
    unsigned char *pixels = new unsigned char[w*h];
    matrix_t<unsigned char> mat_pixels(n);
    // save intensities of pixels into pixels matrix
    image.write(0, 0, w, h, "I", CharPixel, pixels);
    int k = 0;
    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            mat_pixels[i].emplace_back(pixels[k]);
            k++;
        }
    }
    return mat_pixels;
}

template< typename T>
void matrix2image(const matrix_t<T> &matrix, Image &image)
{
    int w = image.columns();
    int h = image.rows();
    // allow image modifying
    image.modifyImage();
    // get values of all pixels
    Quantum *pixels = image.getPixels(0, 0, w, h);
    for (auto &row : matrix) {
        for (auto value : row) {
            *pixels = value * QuantumRange / 255;
            pixels += 2;
        }
    }
    // write changes to pixels
    image.syncPixels();
}

// TODO: void?, a matrix_t<T> &matrix
template< typename T >
void embed(matrix_t<T> &matrix, vector<T> &watermark, int strength, int iteration, int size/*=BLK_SIZE*/)
{
    // threshold values
    double t1 = 3*strength/4;
    double t2 = strength/4;

    for (int i = 0; i < size; ++i) {
        // oposite conditition beacuse 0 is white and 255 is black
        if (watermark.at(i+(size*iteration)) != 0)
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t1;
        else
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t2;
    }
    // return matrix;
}

template< typename T >
vector<T> extract(matrix_t<T> &matrix, int strength, int size/*=BLK_SIZE*/)
{
    vector<T> watermark;
    // threshold values
    double t1 = 3*strength/4;
    double t2 = strength/4;

    // NAOPAK podmienka -> kvoli tomu ze 0-cierna, 255-biela???
    for (int i = 0; i < size; ++i) {
        if (fmod(matrix.at(0).at(i), strength) > ((t1+t2)/2))
            watermark.push_back(255);
        else
            watermark.push_back(0);
    }

    return watermark;
}

template< typename T >
void mat_init(matrix_t<T> &matrix, int size)
{
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            matrix[i].emplace_back(0);
}

// matrix multiplication
// just for n x n matrixes
template< typename T >
matrix_t<T> mat_multiply(matrix_t<T> &a, matrix_t<T> &b, int size/*=BLK_SIZE*/)
{
    matrix_t<T> result(size);
    mat_init(result, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            result.at(i).at(j) = 0;
            for (int k = 0; k < size; ++k)
                result.at(i).at(j) += a.at(i).at(k) * b.at(k).at(j);
        }
    }
    return result;
}

// transpose matrix
template< typename T >
matrix_t<T> mat_transpose(matrix_t<T> &matrix)
{
    matrix_t<T> transposed(matrix.size());
    mat_init(transposed, matrix.size());
    for (int i = 0; i < matrix.size(); ++i)
        for (int j = 0; j < matrix.size(); ++j)
            transposed[j][i] = matrix[i][j];

    return transposed;
}

// compute dot product (inner/scalar product) of two vectors
template< typename T >
double dot_product(const vector<T> &vec1, const vector<T> &vec2)
{
    double sum = 0;
    for (int i = 0; i < vec1.size(); ++i)
        sum += vec1[i] * vec2[i];
    return sum;
}

// normalize vector
// ||x||
template< typename T >
double normalize(const vector<T> &vec)
{
    return sqrt(dot_product(vec, vec));
}

// vector concatenation
template< typename T >
void vector_concat(const vector<T> &vec, vector<T> &result)
{
    result.insert(result.end(), vec.begin(), vec.end());
}

// QR decomposition using Gram-Schmidt process
// https://ocw.mit.edu/courses/mathematics/18-335j-introduction-to-numerical-methods-fall-2010/lecture-notes/MIT18_335JF10_lec10a_hand.pdf
// A = Q*R
// Q is orthogonal matrix
// R is upper triangular matrix
template< typename T >
void gram_schmidt(matrix_t<T> &matrix, matrix3d_t<T> &q_list, matrix3d_t<T> &r_list)
{
    matrix = mat_transpose(matrix);
    matrix_t<double> v(matrix.size());
    matrix_t<double> r(matrix.size());
    matrix_t<double> q(matrix.size());
    mat_init(q, matrix.size());
    mat_init(r, matrix.size());
    v = matrix;
    for (int i = 0; i < matrix.size(); ++i) {
        r[i][i] = normalize(v[i]);
        for (int j = 0; j < matrix.size(); ++j)
            q[i][j] = v[i][j]/r[i][i];
        for (int j = i+1; j < matrix.size(); ++j)
        {
            r[i][j] = dot_product(q[i], v[j]);
            for (int k = 0; k < matrix.size(); ++k)
                v[j][k] -= r[i][j] * q[i][k];
        }
    }

    q = mat_transpose(q);
    matrix = mat_transpose(matrix);
    q_list.emplace_back(q);
    r_list.emplace_back(r);
    // cout << "Q:" << endl;
    // print(q);
    // cout << "R:" << endl;
    // print(r);
    // cout << "--------------test equality----------" << endl;
    // cout << "Original" << endl;
    // print(matrix);
    // cout << "Multiplied" << endl;
    // print(mat_multiply(q, r, matrix.size()));
    // cout << "-------------------------------------" << endl;
}

template< typename T >
void mat_round(matrix_t<T> &matrix)
{
    for (int i = 0; i < matrix.size(); ++i)
        for (int j = 0; j < matrix.size(); ++j)
            matrix[i][j] = round(matrix[i][j]);
}

// peak signal to noise ratio
template< typename T >
double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm)
{
    double psnr = 0;
    double mse = 0; // mean square error

    for (int i = 0; i < orig.size(); ++i)
        for (int j = 0; j < orig.size(); ++j)
            mse += pow((orig[i][j] - wm[i][j]), 2.0);
    mse = mse/(orig.size()*orig.size());
    psnr = 20*log10(255/sqrt(mse));

    return psnr;
}

// normalized correlation
// between original and extracted watermarmark
template< typename T >
double nc(const vector<T> &orig, const vector<T> &extr)
{
    double sum1 = 0;
    double sum2 = 0;

    for (int i = 0; i < orig.size(); ++i) {
        sum1 += (orig[i] * extr[i])/65025;
        sum2 += pow(orig[i]/255, 2.0);
    }

    return sum1/sum2;
}

template< typename T > 
void print(const matrix_t<T> &matrix)
{
    for(const auto& row : matrix) {
        for(double value : row)
            cout << value << " ";
        cout << endl;
    }
    cout << "--------------------" << endl;
}

int main(int argc, char **argv)
{
    InitializeMagick(*argv);

    Image image;
    matrix_t<unsigned char> matrix(image.columns());
    matrix_t<unsigned char> original(image.columns());
    print(matrix);
    // load input image
    image.read(argv[1]);
    matrix = image2matrix(image);
    original = matrix;
    print(matrix);

    // for (int i = 0; i < image.rows(); ++i) {
    //     for (int j = 0; j < image.columns(); ++j) {
    //         matrix.at(i).at(j) = (unsigned char)((i < 4) ? 0 : 255);
    //     }
    // }

    // vector<unsigned char> v = {0,255,0,255,0,255,0,255};
    // embed(matrix, v, atoi(argv[2]), 0, 8);
    // print(matrix);
    // print(original);
    // vector<unsigned char> vv(8);
    // vv = extract(matrix, atoi(argv[2]), 8);
    // for (auto value : vv) {
    //     cout << +value << " ";
    // }
    // cout << endl;


    // print(matrix);
    // cout << "psnr: " << psnr(original, matrix) << endl;
    // cout << "nc: " << nc(v, vv) << endl;

    matrix_t<double> gram_test(original.size());
    for(int i = 0; i < original.size(); ++i)
        for(int j = 0; j < original.size(); ++j)
            gram_test[i].emplace_back(original.at(i).at(j));
    cout << "gram----------" << endl;
    // print(gram_test);
    matrix3d_t<double> q_list;
    matrix3d_t<double> r_list;
    matrix3d_t<double> q_list_out;
    matrix3d_t<double> r_list_out;
    // print(gram_test);
    gram_schmidt(gram_test, q_list, r_list);

    vector<double> vvv = {0,255,0,255,0,255,0,255};
    vector<double> vvvv(8);
    matrix_t<double> gram_embedded(gram_test);
    // q_list[0] = mat_transpose(q_list[0]);
    // cout << "R-list[0]: " << endl;
    // print(r_list[0]);
    // print(q_list[0]);
    // gram_embedded = mat_multiply(q_list[0], r_list[0], 8);
    // print(gram_embedded);
    embed(r_list[0], vvv, atoi(argv[2]), 0, 8);
    // cout << "R-list[0] embedded: " << endl;
    // print(r_list[0]);
    gram_embedded = mat_multiply(q_list[0], r_list[0], 8);

    gram_schmidt(gram_embedded, q_list_out, r_list_out);
    print(gram_test);
    print(gram_embedded);
    mat_round(gram_embedded);
    print(gram_embedded);

    vvvv = extract(r_list_out[0], atoi(argv[2]), 8);
    for (auto value : vvvv) {
        cout << value << " ";
    }
    cout << endl;
    for (auto value : vvv) {
        cout << value << " ";
    }
    cout << endl;

    cout << "psnr: " << psnr(gram_test, gram_embedded) << endl;
    cout << "nc: " << nc(vvv, vvvv) << endl;

    // matrix_t<double> mat{ {12, -51, 4}, {6, 167, -68}, {-4, 24, -41} };
    // // mat_init(mat, 3);
    // // for (int i = 0; i < 3; ++i)
    // // {
    // //     for (int j = 0; j < 3; ++j)
    // //     {
    // //         mat[i][j] = (i+1)*(j+1);
    // //     }
    // // }


    // print(mat);
    // mat = mat_transpose(mat);
    // gram_schmidt(mat, q_list, r_list);
    // print(mat);
    // mat = mat_multiply(mat, mat, 3);
    // print(mat);


    matrix2image(gram_embedded, image);
    image.write("test4.jpg");
    return 0;
}