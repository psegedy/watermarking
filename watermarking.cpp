#include "watermarking.hpp"

using namespace std; 
using namespace Magick;


int main(int argc, char **argv)
{
    InitializeMagick(*argv);

    args_t args = parse_arguments(argc, argv);
    Image image;
    // print(matrix);
    // load input image and watermark
    image.read(get<arg_input>(args));
    

    matrix3d_t<double> blocks(image.columns()*image.rows()/(BLK_SIZE*BLK_SIZE));
    matrix_t<double> matrix(image.columns());
    matrix_t<double> original(image.columns());
    
    matrix3d_t<double> q_list;
    matrix3d_t<double> r_list;
    matrix3d_t<double> q_list_out;
    matrix3d_t<double> r_list_out;
    matrix = image2matrix(image);
    original = matrix;

    // print(matrix);

    // --watermark provided
    if (!get<arg_watermark_in>(args).empty()) {
        Image watermark;
        watermark.read(get<arg_watermark_in>(args));
        matrix_t<double> watermark_m(watermark.columns());
        vector<double> watermark_v;
        vector<double> ext_watermark_v;
        watermark_m = image2matrix(watermark);

        // put watermark to vector
        for (const auto &row : watermark_m)
            for (double value : row)
                watermark_v.emplace_back(value);

        // divide image into blocks
        blocks = get_blocks(matrix);

        // iterate through all blocks
        for (size_t i = 0; i < blocks.size(); ++i) {
            // QR decomposition usinf gram-schmidt process on all blocks
            gram_schmidt(blocks[i], q_list, r_list);
            if (get<arg_embed>(args))
                // embed watermark in R matrix
                embed(r_list[i], watermark_v, get<arg_strength>(args), i);

            if (get<arg_extract>(args))
                // extract watermark from image
                extract(r_list[i], ext_watermark_v, get<arg_strength>(args));

            // multiply Q and R matrixes
            blocks[i] = mat_multiply(q_list[i], r_list[i]);
            // round to 8 bytes
            mat_round(blocks[i]);
        }

        // concatenate blocks
        concat_blocks(matrix, blocks);
        matrix2image(matrix, image);

        // save extracted watermark
        if (get<arg_extract>(args)) {
            // vector to matrix
            for (size_t i = 0, k = 0; i < watermark.columns(); ++i)
                for (size_t j = 0; j < watermark.rows(); ++j, ++k)
                    watermark_m[i][j] = ext_watermark_v[k];

            // save extracted watermark
            matrix2image(watermark_m, watermark);
            watermark.write(get<arg_watermark_out>(args));
        }

        // --attack
        if (get<arg_attack>(args)) {
            // matrix_t<double> embedded = matrix;
            Image attack = image;
            // scale
            Geometry scale75 = Geometry("75x75+0+0%");
            Geometry scale125 = Geometry("125x125+0+0%");

            attack.modifyImage();
            attack.gaussianBlur(1, 1);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "gaussianBlur", get<arg_strength>(args));

            // next attack
            attack = image;
            attack.modifyImage();
            attack.blur();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "blur", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.addNoise(GaussianNoise);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "GaussianNoise", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.addNoise(UniformNoise);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "UniformNoise", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.addNoise(LaplacianNoise);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "LaplacianNoise", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.addNoise(PoissonNoise);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "PoissonNoise", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.medianFilter();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "medianFilter", get<arg_strength>(args));

            // average filter
            attack = image;
            attack.modifyImage();
            const double mask[][3] = {
              {1/9.0, 1/9.0, 1/9.0},
              {1/9.0, 1/9.0, 1/9.0},
              {1/9.0, 1/9.0, 1/9.0}
            };
            attack.convolve(3, &mask[0][0]);
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "averageFilter", get<arg_strength>(args));


            attack = image;
            attack.modifyImage();
            attack.enhance();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "enhance", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.despeckle();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "despeckle", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.sharpen();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "sharpen", get<arg_strength>(args));

            attack = image;
            attack.modifyImage();
            attack.equalize();
            attack.syncPixels();
            do_attack(attack, watermark, watermark_v, "equalize", get<arg_strength>(args));

            Image watermark2 = watermark;
            attack = image;
            attack.modifyImage();
            watermark2.modifyImage();
            watermark2.scale(scale75);
            attack.scale(scale75);
            attack.syncPixels();
            watermark2.syncPixels();
            do_attack(attack, watermark2, watermark_v, "scale75", get<arg_strength>(args), (attack.columns()/(double)image.columns()));

            watermark2 = watermark;
            attack = image;
            attack.modifyImage();
            watermark2.modifyImage();
            watermark2.scale(scale125);
            attack.scale(scale125);
            attack.syncPixels();
            watermark2.syncPixels();
            do_attack(attack, watermark2, watermark_v, "scale125", get<arg_strength>(args), (attack.columns()/(double)image.columns()));
            
            // rotation
            watermark2 = watermark;
            attack = image;
            attack.modifyImage();
            watermark2.modifyImage();
            attack.rotate(5);
            attack.resize(Geometry(560, 560));
            watermark2.resize(Geometry(280, 140));
            attack.syncPixels();
            watermark2.syncPixels();
            do_attack(attack, watermark2, watermark_v, "rotation5", get<arg_strength>(args), (attack.columns()/(double)image.columns()));

            attack = image;
            watermark2 = watermark;
            attack.modifyImage();
            watermark2.modifyImage();
            attack.rotate(10);
            attack.resize(Geometry(600, 600));
            watermark2.resize(Geometry(300, 150));
            attack.syncPixels();
            watermark2.syncPixels();
            do_attack(attack, watermark2, watermark_v, "rotation10", get<arg_strength>(args), (attack.columns()/(double)image.columns()));

            // JPEG
            attack = image;
            attack.quality(70);
            attack.write("attacks/jpeg70.jpg");
            attack.read("attacks/jpeg70.jpg");
            do_attack(attack, watermark, watermark_v, "jpeg70", get<arg_strength>(args));

            attack = image;
            attack.quality(50);
            attack.write("attacks/jpeg50.jpg");
            attack.read("attacks/jpeg50.jpg");
            do_attack(attack, watermark, watermark_v, "jpeg50", get<arg_strength>(args));

            attack = image;
            attack.quality(30);
            attack.write("attacks/jpeg30.jpg");
            attack.read("attacks/jpeg30.jpg");
            do_attack(attack, watermark, watermark_v, "jpeg30", get<arg_strength>(args));

            attack = image;
            attack.quality(20);
            attack.write("attacks/jpeg20.jpg");
            attack.read("attacks/jpeg20.jpg");
            do_attack(attack, watermark, watermark_v, "jpeg20", get<arg_strength>(args));

        }

        // --statistics
        if (get<arg_stat>(args)) {
            // PSNR and NC
            cout << "PSNR: " << psnr(original, matrix) << endl;
            cout << "NC: " << nc(watermark_v, ext_watermark_v) << endl;
        }
    }
    
    // set compression level
    image.quality(get<arg_quality>(args));
    // save watermarked image
    image.write(get<arg_output>(args));

    return 0;
}

// FUNCTIONS

void print_help()
{
    cout << "--help - show this help" << endl;
    cout << "--input FILE - define input image" << endl;
    cout << "--output FILE - define output file - default: watermarked.jpg" << endl;
    cout << "--wm-input FILE - watermark in image file" << endl;
    cout << "--wm-output FILE - file where save extracted watermark - default: extracted.jpg" << endl;
    cout << "--strength NUMBER - set watermarking strength" << endl;
    cout << "--quality NUMBER - set compression level for output (watermarked) image in range <0-100>" << endl;
    cout << "--embed - embedding" << endl;
    cout << "--extract - extraction" << endl;
    cout << "--attack - do all attacks, save images to attacks/ directory" << endl;
    cout << "--statistics - show NC of extracted and original watermark and PSNR of original and watermarked image, use with --extract or without --embed and --extract" << endl;
    cout << endl;
    cout << "--input and --wm-input are required" << endl;
    cout << "If --embed and --extract are not specified, it runs embedding of watermark and then extraction of embedded watermark" << endl;
}

void do_attack(Image &attack, Image &watermark, vector<double> &watermark_v, string name, int strength, double scale/*=1*/)
{
    matrix_t<double> matrix_attack = image2matrix(attack);
    matrix3d_t<double> q_list;
    matrix3d_t<double> r_list;
    matrix3d_t<double> blocks(attack.columns()*attack.rows()/(BLK_SIZE*BLK_SIZE));
    vector<double> ext_watermark_v;
    matrix_t<double> watermark_m;
    watermark_m = image2matrix(watermark);

    // divide image into blocks
    blocks = get_blocks(matrix_attack);
    for (size_t i = 0; i < blocks.size(); ++i) {
        // QR decomposition usinf gram-schmidt process on all blocks
        gram_schmidt(blocks[i], q_list, r_list);
        // extract watermark from image
        extract(r_list[i], ext_watermark_v, strength);
    }
    // vector to matrix
    for (size_t i = 0, k = 0; i < watermark.columns(); ++i)
        for (size_t j = 0; j < watermark.rows(); ++j, ++k)
            watermark_m[i][j] = ext_watermark_v[k];

    // save extracted watermark
    matrix2image(watermark_m, watermark);
    watermark.write("attacks/" + name + "-wm.jpg");
    attack.write("attacks/" + name + ".jpg");
    cout << name << endl;

    // resize extracted watermark to default watermark size to compare original and extracted
    if (scale != 1) {
        watermark.modifyImage();
        watermark.resize(Geometry(watermark.columns()/scale, watermark.rows()/scale));
        watermark.syncPixels();

        watermark_m = image2matrix(watermark);
        // put watermark to vector
        ext_watermark_v.clear();
        for (const auto &row : watermark_m)
            for (double value : row)
                ext_watermark_v.emplace_back(value);
    }
    cout << "NC: " << nc(watermark_v, ext_watermark_v) << endl;
    cout << "-----------------" << endl;
}

// parse command line arguments
args_t parse_arguments(int argc, char **argv)
{
    args_t args;
    // set some defaults
    get<arg_output>(args) = OUTPUT;
    get<arg_watermark_out>(args) = WM_OUTPUT;
    get<arg_quality>(args) = QUALITY;
    get<arg_strength>(args) = STRENGTH;
    get<arg_embed>(args) = false;
    get<arg_extract>(args) = false;
    get<arg_attack>(args) = false;
    get<arg_stat>(args) = false;

    if (argc == 1) {
        print_help();
        exit(2);
    }
    // vector of argument regexes 
    vector<regex> regexes;
    for (string name : arg_names) {
        regexes.emplace_back(name);
    }

    for (int i = 1; i < argc; ++i) {
        // --help
        if (regex_match(argv[i], regexes[0])) {
            print_help();
            exit(0);
        }
        // other arguments
        for (size_t j = 1; j < regexes.size(); ++j) {
            if (regex_match(argv[i], regexes[j])) {
                if (i < arg_quality && argc < (i+1)) {
                    print_help();
                    exit(2);
                }
                // write value to tuple
                switch (j) {
                    case arg_input :
                        get<arg_input>(args) = argv[i+1];
                        break;
                    case arg_output :
                        get<arg_output>(args) = argv[i+1];
                        break;
                    case arg_watermark_in :
                        get<arg_watermark_in>(args) = argv[i+1];
                        break;
                    case arg_watermark_out :
                        get<arg_watermark_out>(args) = argv[i+1];
                    case arg_strength :
                        get<arg_strength>(args) = atoi(argv[i+1]);
                        break;
                    case arg_quality :
                        get<arg_quality>(args) = atoi(argv[i+1]);
                        break;
                    case arg_embed :
                        get<arg_embed>(args) = true;
                        break;
                    case arg_extract :
                        get<arg_extract>(args) = true;
                        break;
                    case arg_attack :
                        get<arg_attack>(args) = true;
                        break;
                    case arg_stat :
                        get<arg_stat>(args) = true;
                        break;
                    default :
                        break;
                }
            }
        }
    }
    if (get<arg_input>(args).empty()) {
        cerr << "--input <filename> is required" << endl << "use --help for more information" << endl;
        exit(2);
    }
    if (get<arg_stat>(args) == true || get<arg_embed>(args) == true || get<arg_extract>(args) == true) {
        if (get<arg_watermark_in>(args).empty()) {
            cerr << "--wm-input is required with --statistics, --embed or --extract" << endl << "use --help for more information" << endl;
            exit(2);
        }
    }
    if (get<arg_embed>(args) == false && get<arg_extract>(args) == false) {
        get<arg_embed>(args) = true;
        get<arg_extract>(args) = true;
    }
    if (get<arg_attack>(args) == true)
        if (get<arg_embed>(args) == false || get<arg_extract>(args) == false) {
            cerr << "There should be --embed and --extract when using --attack" << endl;
            exit(2);
        }

    return args;
}

// get BLK_SIZE x BLK_SIZE blocks from matrix
template< typename T >
matrix3d_t<T> get_blocks(const matrix_t<T> &matrix, int n/*=BLK_SIZE*/)
{
    matrix_t<T> block(n);
    matrix3d_t<T> blocks;
    for (size_t i = 0; i < matrix.size(); i += n)
        for (size_t j = 0; j < matrix.size(); j += n) {
            block.clear();
            block.resize(n);
            for(int k = 0; k < n; ++k)
                for(int l = 0; l < n; ++l)
                    block[k].emplace_back(matrix.at(k+i).at(l+j));
            blocks.emplace_back(block);
        }

    return blocks;
}

// Concatenate blocks
template< typename T >
void concat_blocks(matrix_t<T> &result, matrix3d_t<T> &blocks)
{
    size_t i = 0;
    size_t j = 0;
    for(const auto block : blocks) {
        for (size_t k = 0; k < block.size(); ++k)
            for (size_t l = 0; l < block.size(); ++l)
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

// transform Image object into matrix_t type
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

// transform matrix of image pixels to Image
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

// Watermark embedding algorithm
template< typename T >
void embed(matrix_t<T> &matrix, vector<T> &watermark, int strength, int iteration, int size/*=BLK_SIZE*/)
{
    // threshold values
    double t1 = 3*strength/4.0;
    double t2 = strength/4.0;

    for (int i = 0; i < size; ++i) {
        // inverse conditition beacuse 0 is white and 255 is black
        if (watermark.at(i+(size*iteration)) != 0)
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t1;
        else
            matrix.at(0).at(i) = matrix.at(0).at(i) - fmod(matrix.at(0).at(i), strength)  + t2;
    }
}

// Algorithm for watermark extraction
template< typename T >
void extract(matrix_t<T> &matrix, vector<T> &watermark, int strength, int size/*=BLK_SIZE*/)
{
    // threshold values
    double t1 = 3*strength/4.0;
    double t2 = strength/4.0;

    for (int i = 0; i < size; ++i) {
        if (fmod(matrix.at(0).at(i), strength) > ((t1+t2)/2.0))
            watermark.push_back(255);
        else
            watermark.push_back(0);
    }
}

// Initialize matrix with zeroes
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
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            transposed[j][i] = matrix[i][j];

    return transposed;
}

// compute dot product (inner/scalar product) of two vectors
template< typename T >
double dot_product(const vector<T> &vec1, const vector<T> &vec2)
{
    double sum = 0;
    for (size_t i = 0; i < vec1.size(); ++i)
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

    for (size_t i = 0; i < matrix.size(); ++i) {
        r[i][i] = normalize(v[i]);  // normalize vector
        for (size_t j = 0; j < matrix.size(); ++j) {
            if (r[i][i] == 0)   // avoid division by zero
                q[i][j] = 0;
            else
                q[i][j] = v[i][j]/r[i][i];
        }
        for (size_t j = i+1; j < matrix.size(); ++j) {
            r[i][j] = dot_product(q[i], v[j]);
            for (size_t k = 0; k < matrix.size(); ++k)
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

// round matrix values
template< typename T >
void mat_round(matrix_t<T> &matrix)
{
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            matrix[i][j] = round(matrix[i][j]);
}

// peak signal to noise ratio
template< typename T >
long double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm)
{
    double psnr = 0;
    double mse = 0; // mean square error

    for (size_t i = 0; i < orig.size(); ++i)
        for (size_t j = 0; j < orig.size(); ++j)
            mse += pow((orig[i][j] - wm[i][j]), 2.0);

    mse = mse/(orig.size()*orig.size());
    psnr = 20*log10(255/sqrt(mse));

    return psnr;
}

// normalized correlation
// between original and extracted watermarmark
template< typename T >
double nc(vector<T> &orig, vector<T> &extr)
{
    double sum1 = 0;
    double sum2 = 0;

    for (size_t i = 0; i < orig.size(); ++i) {
        if (orig[i] == 0)
            orig[i] = -1;
        else
            orig[i] = 1;
        if (extr[i] == 0)
            extr[i] = -1;
        else
            extr[i] = 1;
            
        sum1 += (orig[i] * extr[i]);
        sum2 += pow(orig[i], 2.0);
    }

    return sum1/sum2;
}

// print matrix
// for debugging purposes
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