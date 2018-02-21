# Watermarking as Document Protection
Bachelor thesis by Patrik Segedy and supervised by Doc. Dr. Ing. Petr Hanáček at Brno University of Technology, Faculty of Information Technology, Department of Intelligent Systems

This thesis is dealing with document protection using a digital watermarking. First, a watermark characteristics are presented. Then, different usages of the watermarking are discussed. Next part of the thesis is dedicated to current development in watermarking field. It is aimed at various principles of watermark embedding into different multimedia types. Subsequently, a watermarking scheme for still images is proposed and implemented. Finally, the watermarking scheme undergoes attacks, which should damage embedded watermark. Attacked watermark is then extracted and compared to the embedded watermark.

Implemented algorithm was presented in an article[1], where authors are using [QR Decomposition](https://en.wikipedia.org/wiki/QR_decomposition) for embedding and extracting the watermark to/from greyscale still image.

## How to compile
Program is using [ImageMagick](http://ww.imagemagick.org) library. Thanks to this library, it is possible to use over 200 image formats.
You will need *ImageMagick* in version *> 7.0.5-4* 

To check installed vesion of ImageMagick use:
```
identify -version
```

Before compiling it is neccessary to set environmental variable `PKG_CONFIG_PATH`:
``` bash
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
```

Maybe you will need to configure dynamic linker run-time bindings
``` bash
sudo ldconfig /usr/local/lib
```

Next step is simply run Makefile:
```
make
```
which create executable file `watermarking` and `attacks` folder, where will be stored images after attacks (run atacks with `--attack` option, see options below)

Run the application as:
``` bash
./watermarking --input lena.jpg --wm-input fit.pbm
```

## Options

* --help - show this help
* --input FILE - define input image
* --output FILE - define output file - default: watermarked.jpg
* --wm-input FILE - watermark in image file
* --wm-output FILE - file where save extracted watermark - default: extracted.jpg
* --strength NUMBER - set watermarking strength
* --quality NUMBER - set compression level for output (watermarked) image in range <0-100>
* --embed - embedding
* --extract - extraction
* --attack - do all attacks, save images to attacks/ directory
* --statistics - show NC of extracted and original watermark and PSNR of original and watermarked image, use with --extract or without --embed and --extract

* --input and --wm-input are required
If --embed and --extract are not specified, it runs embedding of watermark and then extraction of embedded watermark

## References
The tesis is available at https://dspace.vutbr.cz/xmlui/handle/11012/69833 (in Slovak language)

[1] NADERAHMADIAN, Yashar a Saied HOSSEINI-KHAYAT. Fast and robust watermarking in still images based on QR decomposition. Multimedia Tools and Applications [online]. 2014, 72(3), 2597-2618 [cit. 2017-05-01]. DOI: 10.1007/s11042-013-1559-9. ISSN 13807501. Available from: http://link.springer.com/10.1007/s11042-013-1559-9
