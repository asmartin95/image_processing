#include <stdio.h>                          /* Sobel.c */
#include <stdlib.h> 
#include <math.h>

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
int masky[3][3] = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
double ival[256][256], maxival;

main(argc, argv)
int argc;
char **argv;
{
	int i, j, p, q, mr, sum1, sum2;
	double threshold;
	FILE *fo1, *fo2, *fo3, *fp1, *fopen();
	char *foobar;

	//input file
	argc--; argv++;
	foobar = *argv;
	fp1 = fopen(foobar, "rb");

	//threshold input
	argc--; argv++;
	foobar = *argv;
	threshold = atof(foobar);

	//Magnitude output file
	fo1 = fopen("sobelMagnitude.pgm", "wb");

	//HiThreshold output file
	fo2 = fopen("sobelHiThreshold.pgm", "wb");
	
	//LoThreshold output file
	fo3 = fopen("sobelLoThreshold.pgm", "wb");

	//Read in header info from input pic (to prevent wrapping)
	int rowIn, colIn, maxValIn;
	fscanf(fp1, "%s %d %d %d", foobar, &rowIn, &colIn, &maxValIn);

	//PGM headers
	fprintf(fo1, "P5\t");
	fprintf(fo1, "%d %d\t", rowIn, colIn);
	fprintf(fo1, "255\n");
	fprintf(fo2, "P5\n");
	fprintf(fo2, "%d %d\n", rowIn, colIn);
	fprintf(fo2, "255\n");
	fprintf(fo3, "P5\n");
	fprintf(fo3, "%d %d\n", rowIn, colIn);
	fprintf(fo3, "255\n");

	//Make Hi Threshold 140% of entered threshold and Lo Threshold 60%
	double HiThreshold = threshold * 1.3;
	double LoThreshold = threshold * 0.7;

	//build picture array
	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < 256; j++)
		{
			pic[i][j] = getc(fp1);
			pic[i][j] &= 0377;
		}
	}

	//apply the x and y mask and store it in the outpic array
	mr = 1;
	for (i = mr; i < 256 - mr; i++)
	{
		for (j = mr; j < 256 - mr; j++)
		{
			sum1 = 0;
			sum2 = 0;
			for (p = -mr; p <= mr; p++)
			{
				for (q = -mr; q <= mr; q++)
				{
					sum1 += pic[i + p][j + q] * maskx[p + mr][q + mr];
					sum2 += pic[i + p][j + q] * masky[p + mr][q + mr];
				}
			}
			outpicx[i][j] = sum1;
			outpicy[i][j] = sum2;
		}
	}

	//calculate magnitude and store in ival array
	maxival = 0;
	for (i = mr; i < 256 - mr; i++)
	{
		for (j = mr; j<256 - mr; j++)
		{
			ival[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) +
				(outpicy[i][j] * outpicy[i][j])));
			if (ival[i][j] > maxival)
				maxival = ival[i][j];

		}
	}

	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < 256; j++)
		{
			ival[i][j] = (ival[i][j] / maxival) * 255;

			//magnitude
			fprintf(fo1, "%c", (char)((int)(ival[i][j])));
			
			//HiThreshold
			if (ival[i][j] > HiThreshold) {
				fprintf(fo2, "%c", (char)((int)(255)));
			}
			else if (ival[i][j] < HiThreshold){
				fprintf(fo2, "%c", (char)((int)(0)));
			}

			//LoThreshold
			if (ival[i][j] > LoThreshold) {
				fprintf(fo3, "%c", (char)((int)(255)));
			}
			else if (ival[i][j] < LoThreshold){
				fprintf(fo3, "%c", (char)((int)(0)));
			}

		}
	}

}