#include <vector>
#include <Magick++.h> 
#include <iostream>
#include <cstdlib> //atoi
#include <cmath>
#include <regex>
#include <string>
#include <tuple>

#define BLK_SIZE 8
#define OUTPUT "watermarked.jpg"
#define WM_OUTPUT "extracted.jpg"
#define QUALITY 100
#define STRENGTH 24

using namespace std; 
using namespace Magick;

// matrix_t<T> type
template< typename T > using matrix_t = vector< vector<T> >;

// matrix3d_t<T> type
template< typename T > using matrix3d_t = vector< matrix_t<T> >;

enum arg_enum { arg_help = 0, arg_input, arg_output, arg_watermark_in, arg_watermark_out, arg_strength, arg_quality, arg_embed, arg_extract, arg_attack, arg_stat };
vector<string> arg_names = { "--help", "--input", "--output", "--wm-input", "--wm-output", "--strength", "--quality", "--embed", "--extract", "--attack", "--statistics" };
// args_t type for command line arguments
using args_t = tuple<bool, string, string, string, string, int, int, bool, bool, bool, bool>;

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
void mat_round(matrix_t<T> &matrix);
template< typename T >
matrix_t<T> mat_multiply(matrix_t<T> &a, matrix_t<T> &b, int size=BLK_SIZE);
template< typename T >
matrix_t<T> mat_transpose(matrix_t<T> &matrix);
template< typename T >
double dot_product(const vector<T> &vec1, const vector<T> &vec2);
template< typename T >
double normalize(const vector<T> &vec);
template< typename T >
void gram_schmidt(matrix_t<T> &matrix, matrix3d_t<T> &q_list, matrix3d_t<T> &r_list);
template< typename T >
long double psnr(const matrix_t<T> &orig, const matrix_t<T> &wm);
template< typename T >
double nc(vector<T> &orig, vector<T> &extr);
template< typename T > 
void print(const matrix_t<T> &matrix);
args_t parse_arguments(int argc, char **argv);
void do_attack(Image &attack, Image &watermark, vector<double> &watermark_v, string name, int strength, double scale=1);
void print_help();
