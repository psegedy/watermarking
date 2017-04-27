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
matrix3d_t<T> get_blocks(const matrix_t<T> &matrix, int n=BLK_SIZE);
template< typename T >
void concat_blocks(matrix_t<T> &result, matrix3d_t<T> &blocks);
matrix_t<double> image2matrix(Image &image);
template< typename T>
void matrix2image(const matrix_t<T> &matrix, Image &image);
template< typename T >
void embed(matrix_t<T> &matrix, vector<T> &watermark, int strength, int iteration, int size=BLK_SIZE);
template< typename T >
void extract(matrix_t<T> &matrix, vector<T> &watermark, int strength, int size=BLK_SIZE);
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
long double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm);
template< typename T >
double nc(const vector<T> &orig, const vector<T> &extr);
template< typename T > 
void print(const matrix_t<T> &matrix);


// get BLK_SIZE x BLK_SIZE blocks from matrix
template< typename T >
matrix3d_t<T> get_blocks(const matrix_t<T> &matrix, int n/*=BLK_SIZE*/)
{
    matrix_t<T> block(n);
    matrix3d_t<T> blocks;
    for (int i = 0; i < matrix.size(); i += n)
        for (int j = 0; j < matrix.size(); j += n) {
            block.clear();
            block.resize(n);
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
        for (int k = 0; k < block.size(); ++k)
            for (int l = 0; l < block.size(); ++l)
                result.at(k+i).at(l+j) = block.at(k).at(l);
        if (j < result.size() - BLK_SIZE) {
            j += BLK_SIZE;
        }
        else {
            j = 0;
            if (i < result.size() - BLK_SIZE)
                i += BLK_SIZE;
        }
    }
}

// transform Image objet into matrix_t type
// every value is in range <0-255>
//   0     - black
//   255   - white
matrix_t<double> image2matrix(Image &image)
{
    int w = image.columns();
    int h = image.rows();
    int n = w; // image is n x n
    unsigned char *pixels = new unsigned char[w*h];
    matrix_t<double> mat_pixels(n);
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
            pixels += image.channels();
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
        // inverse conditition beacuse 0 is white and 255 is black
        if (watermark.at(i+(size*iteration)) != 0)
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t1;
        else
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t2;
    }
}

template< typename T >
void extract(matrix_t<T> &matrix, vector<T> &watermark, int strength, int size/*=BLK_SIZE*/)
{
    // threshold values
    double t1 = 3*strength/4;
    double t2 = strength/4;

    for (int i = 0; i < size; ++i) {
        if (fmod(matrix.at(0).at(i), strength) > ((t1+t2)/2.0))
            watermark.push_back(255);
        else
            watermark.push_back(0);
    }
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
    // transpose input matrix
    matrix = mat_transpose(matrix);
    matrix_t<double> v(matrix.size());
    matrix_t<double> r(matrix.size());
    matrix_t<double> q(matrix.size());
    // initialize matrixes to zeroes
    mat_init(q, matrix.size());
    mat_init(r, matrix.size());
    v = matrix;

    for (int i = 0; i < matrix.size(); ++i) {
        r[i][i] = normalize(v[i]);  // normalize vector
        for (int j = 0; j < matrix.size(); ++j) {
            if (r[i][i] == 0)   // avoid division by zero
                q[i][j] = 0;
            else
                q[i][j] = v[i][j]/r[i][i];
        }
        for (int j = i+1; j < matrix.size(); ++j)
        {
            r[i][j] = dot_product(q[i], v[j]);
            for (int k = 0; k < matrix.size(); ++k)
                v[j][k] -= r[i][j] * q[i][k];
        }
    }

    // transpose Q matrix
    q = mat_transpose(q);
    // transpose input matrix to original state
    matrix = mat_transpose(matrix);
    // store computed R,Q matrixes for all blocks
    q_list.emplace_back(q);
    r_list.emplace_back(r);
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
long double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm)
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
    for(const auto &row : matrix) {
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
    Image watermark;
    // print(matrix);
    // load input image and watermark
    image.read(argv[1]);
    watermark.read(argv[2]);

    matrix3d_t<double> blocks(image.columns()*image.rows()/(BLK_SIZE*BLK_SIZE));
    matrix_t<double> matrix(image.columns());
    matrix_t<double> original(image.columns());
    matrix_t<double> watermark_m(watermark.columns());
    vector<double> watermark_v;
    vector<double> ext_watermark_v;
    matrix3d_t<double> q_list;
    matrix3d_t<double> r_list;
    matrix3d_t<double> q_list_out;
    matrix3d_t<double> r_list_out;
    matrix = image2matrix(image);
    watermark_m = image2matrix(watermark);
    original = matrix;

    // put watermark to vector
    for (const auto &row : watermark_m)
        for (double value : row)
            watermark_v.emplace_back(value);

    // print(matrix);

    // divide image into blocks
    blocks = get_blocks(matrix);

    // iterate through all blocks
    for (int i = 0; i < blocks.size(); ++i) {
        // QR decomposition usinf gram-schmidt process on all blocks
        gram_schmidt(blocks[i], q_list, r_list);
        
        // embed watermark in R matrix
        embed(r_list[i], watermark_v, atoi(argv[3]), i);

        // extract watermark from image
        extract(r_list[i], ext_watermark_v, atoi(argv[3]));

        // multiply Q and R matrixes
        blocks[i] = mat_multiply(q_list[i], r_list[i]);
        // round to 8 bytes
        mat_round(blocks[i]);
    }
    

    // concatenate blocks
    concat_blocks(matrix, blocks);
    matrix2image(matrix, image);
    // image.quality(30);
    // save watermarked image
    image.write("watermarked.jpg");


    // extract watermark from image
    // for (auto r : r_list)
    //     extract(r, ext_watermark_v, atoi(argv[3]));

    // PSNR and NC
    cout << "psnr: " << psnr(original, matrix) << endl;
    cout << "nc: " << nc(watermark_v, ext_watermark_v) << endl;
    int count = 0;
    for (int i = 0; i < watermark_v.size(); ++i)
        if (watermark_v[i] != ext_watermark_v[i])
            count++;
    cout << "count: " << count << endl;

    // vector to matrix
    for (int i = 0, k = 0; i < watermark.columns(); ++i)
        for (int j = 0; j < watermark.rows(); ++j, ++k)
            watermark_m[i][j] = ext_watermark_v[k];

    // save extracted watermark
    matrix2image(watermark_m, watermark);
    watermark.write("ext_watermark.jpg");

    return 0;
}