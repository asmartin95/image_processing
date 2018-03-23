#include <stdio.h>                          /* Canny.c */
#include <stdlib.h> 
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

int    pic[PICSIZE][PICSIZE];
int histogram[PICSIZE] = { 0 };
int peakFlag[PICSIZE][PICSIZE] = { 0 };
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE];
double convy[PICSIZE][PICSIZE];
double peak[PICSIZE][PICSIZE];
double fin[PICSIZE][PICSIZE];
double mag[PICSIZE][PICSIZE];

main(argc, argv)
int argc;
char **argv;
{
	int     i, j, x, y, p, q, mr, centx, centy;
	double  maskval, sum1, sum2, sig, maxival, HI, LO;
	FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
	char    *foobar;

	//input file
	argc--; argv++;
	foobar = *argv;
	fp1 = fopen(foobar, "rb");

	//Read in header info from input pic (to prevent wrapping)
	int rowIn, colIn, maxValIn;
	fscanf(fp1, "%s %d %d %d", foobar, &rowIn, &colIn, &maxValIn);

	//set sigma
	argc--; argv++;
	foobar = *argv;
	sig = atof(foobar);

	//Magnitude image
	fo1 = fopen("magnitude.pgm", "wb");
	fprintf(fo1, "P5\t");
	fprintf(fo1, "%d %d\t", rowIn, colIn);
	fprintf(fo1, "255\n");

	//Peaks image
	fo2 = fopen("peaks.pgm", "wb");
	fprintf(fo2, "P5\n");
	fprintf(fo2, "%d %d\n", rowIn, colIn);
	fprintf(fo2, "255\n");

	//Final Edges image
	fo3 = fopen("finalEdge.pgm", "wb");
	fprintf(fo3, "P5\n");
	fprintf(fo3, "%d %d\n", rowIn, colIn);
	fprintf(fo3, "255\n");

	mr = (int)(sig * 3);
	centx = (MAXMASK / 2);
	centy = (MAXMASK / 2);

	//build picture array
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			pic[i][j] = getc(fp1);
		}
	}

	//build mask
	for (p = -mr; p <= mr; p++) {
		for (q = -mr; q <= mr; q++) {
			//Gaussian filter formula
			maskval = p * exp(-((p*p + q*q) / (2 * (sig * sig))));
			(maskx[p + centy][q + centx]) = maskval;
			maskval = q * exp(-((q*q + p*p) / (2 * (sig * sig))));
			(masky[p + centy][q + centx]) = maskval;
		}
	}

	//do convolution
	for (i = mr; i <= 256 - mr; i++) {
		for (j = mr; j <= 256 - mr; j++) {
			sum1 = 0;
			sum2 = 0;
			for (p = -mr; p <= mr; p++) {
				for (q = -mr; q <= mr; q++) {
					sum1 += pic[i + p][j + q] * maskx[p + centy][q + centx];
					sum2 += pic[i + p][j + q] * masky[p + centy][q + centx];
				}
			}
			outpicx[i][j] = sum1;
			outpicy[i][j] = sum2;
			convx[i][j] = sum1;
			convy[i][j] = sum2;
		}
	}

	//calculate magnitude and store in mag array
	maxival = 0;
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			mag[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) + (outpicy[i][j] * outpicy[i][j])));
			if (mag[i][j] > maxival) maxival = mag[i][j];
		}
	}

	//print magnitude
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			mag[i][j] = (mag[i][j] / maxival) * 255;
			fprintf(fo1, "%c", (char)((int)(mag[i][j])));
		}
	}

	//Find peak map
	double slope = 0.0;
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			if ((convx[i][j]) == 0.0) {
				convx[i][j] = .0000001;
			}

			slope = convy[i][j] / convx[i][j];

			//Horizontal slope (check positions N/S)
			if ((slope <= tan(22.5 * M_PI / 180)) && (slope > tan(-22.5 * M_PI / 180))) {
				if ((mag[i][j] > mag[i - 1][j]) && (mag[i][j] > mag[i + 1][j])) {
					peak[i][j] = 255;
					peakFlag[i][j] = 1;
				}
			}

			//SW/NE slope (check positions SE + NW)
			else if ((slope <= tan(67.5 * M_PI / 180)) && (slope > tan(22.5 * M_PI / 180))) {
				if ((mag[i][j] > mag[i + 1][j + 1]) && (mag[i][j] > mag[i - 1][j - 1])) {
					peak[i][j] = 255;
					peakFlag[i][j] = 1;
				}
			}

			//SE/NW test (check positions SW + NE)
			else if ((slope <= tan(-22.5 * M_PI / 180)) && (slope > tan(-67.5 * M_PI / 180))) {
				if ((mag[i][j] > mag[i - 1][j + 1]) && (mag[i][j] > mag[i + 1][j - 1])) {
					peak[i][j] = 255;
					peakFlag[i][j] = 1;
				}
			}

			//Vertical test (check positions E + W)
			else {
				if ((mag[i][j] > mag[i][j - 1]) && (mag[i][j] > mag[i][j + 1])) {
					peak[i][j] = 255;
					peakFlag[i][j] = 1;
				}
			}
		}
	}

	//Output peak image
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			fprintf(fo2, "%c", (char)((int)(peak[i][j])));
		}
	}

	//Build Histogram
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			(histogram[(int)(mag[i][j])])++;
		}
	}

	//calculate HI/LO
	double percent = 0.05;
	double cutoff = percent * rowIn * colIn;
	double areaOfTops = 0.0;
	for (i = 255; i >= 0; i--) {
		areaOfTops += histogram[i];
		if (areaOfTops > cutoff) {
			HI = i;
			LO = 0.35*HI;
			break;
		}
	}

	//Build final image with mag and peakflag arrays
	int moreToDo;
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			//inefficient implementation
			if (peakFlag[i][j] == 1) {
				//HI
				if (mag[i][j] > HI) {
					peakFlag[i][j] = 0;
					fin[i][j] = 255;
				}
				//LO
				else if (mag[i][j] < LO) {
					peakFlag[i][j] = 0;
					fin[i][j] = 0;
				}
				//Middle loop
				//inefficient but easy to implement method
				else {
					moreToDo = 1;
					while (moreToDo == 1) {
						moreToDo = 0;
						for (x = i; x < 256; x++) {
							for (y = j; y < 256; y++) {
								if (peakFlag[x][y] == 1) {
									for (p = -1; p <= 1; p++) {
										for (q = -1; q <= 1; q++) {
											if (fin[x + p][y + q] == 255) {
												peakFlag[x][y] = 0;
												fin[x][y] = 255;
												moreToDo = 1;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	//Output final image
	for (i = 0; i < 256; i++) {
		for (j = 0; j < 256; j++) {
			fprintf(fo3, "%c", (char)((int)(fin[i][j])));
		}
	}
}

