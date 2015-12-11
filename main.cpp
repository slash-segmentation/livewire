/**

Livewire Qt Demo - Demo for running the Livewire algorithm
Copyright (C) 2011  Jeffrey Bush  jeff@coderforlife.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

**/

#include "QtPanel.h"

#include <QtGui/QApplication>
#include <QImage>

#define DEFAULT_IMAGE "test.png" // or "test-color-0.jpg" for a color image

//#define RUN_WEIGHT_TESTS

#ifdef RUN_WEIGHT_TESTS
void run_all_weights();
#endif

int main(int argc, char *argv[])
{
#ifdef RUN_WEIGHT_TESTS
	run_all_weights();
	return 0;
#else
	QApplication app(argc, argv);
	Livewire::QtPanel w;
	QImage img((argc > 1) ? argv[1] : DEFAULT_IMAGE);
	w.SetImage(img);
	w.show();
	return app.exec();
#endif
}

#ifdef RUN_WEIGHT_TESTS
#include "Weights.h"
void run_all_weights()
{
	QImage img(DEFAULT_IMAGE);

	const uint w = img.width(), h = img.height();
	const uint cx = w / 2, cy = h / 2, min_room = w > h ? w - cx : h - cy; // given to CalculateRegion to force the whole image to calculate
	const bool gray = img.format() == QImage::Format_Indexed8 && img.isGrayscale();
	img = gray ? img.copy() : img.convertToFormat(QImage::Format_RGB32);

	Livewire::Weights weights;
	weights.SetImage(img.bits(), w, h, gray ? Livewire::Weights::GrayscaleByte : Livewire::Weights::RGB, img.bytesPerLine());

	static const Livewire::Weights::PixelReductionMethod px_reds[] = { Livewire::Weights::NoPixelReduction,
		Livewire::Weights::Median2pxWindow, Livewire::Weights::Median3pxWindow, Livewire::Weights::Median4pxWindow, Livewire::Weights::Median5pxWindow,
		Livewire::Weights::Mean2pxWindow, Livewire::Weights::Mean3pxWindow, Livewire::Weights::Mean4pxWindow, Livewire::Weights::Mean5pxWindow,
		Livewire::Weights::Gaussian3pxWindow, Livewire::Weights::Gaussian4pxWindow, Livewire::Weights::Gaussian5pxWindow,
	};
	static const char* px_red_strs[] = { "none", "med2", "med3", "med4", "med5", "mean2", "mean3", "mean4", "mean5", "gauss3", "gauss4", "gauss5", };
	static const int n_pr = sizeof(px_reds) / sizeof(px_reds[0]);

	Livewire::Weights::NoiseReductionMethod noise_reds[] = { Livewire::Weights::NoNoiseReduction,
		Livewire::Weights::MedianFilter3pxWindow, Livewire::Weights::MedianFilter5pxWindow,
		Livewire::Weights::MeanFilter3pxWindow, Livewire::Weights::MeanFilter5pxWindow,
		Livewire::Weights::GaussianFilter3pxWindow, Livewire::Weights::GaussianFilter5pxWindow,
	};
	static const char* noise_red_strs[] = { "none", "med3", "med5", "mean3", "mean5", "gauss3", "gauss5", };
	static const int n_nr = sizeof(noise_reds) / sizeof(noise_reds[0]);

	Livewire::Weights::AccentuationMethod accentuates[] = { Livewire::Weights::NoAccentuation, Livewire::Weights::Sigmoid, };
	static const char* accentuate_strs[] = { "none", "sigmoid", };
	static const int n_a = sizeof(accentuates) / sizeof(accentuates[0]);

	Livewire::Weights::EdgeDetectionMethod edge_detects[] = { Livewire::Weights::NoEdgeDetection,
		Livewire::Weights::Sobel3, Livewire::Weights::Sobel5, Livewire::Weights::Scharr,
		Livewire::Weights::Canny,
	};
	static const char* edge_detect_strs[] = { "none", "sobel3", "sobel5", "scharr", "canny", };
	static const int n_ed = sizeof(edge_detects) / sizeof(edge_detects[0]);

	Livewire::Weights::Settings settings;
	char name[1024];
	for (int pr = 0; pr < n_pr; ++pr)
	{
		settings.PixelReduction = px_reds[pr];
		for (int nr = 0; nr < n_nr; ++nr)
		{
			settings.NoiseReduction = noise_reds[nr];
			for (int a = 0; a < n_a; ++a)
			{
				settings.Accentuation = accentuates[a];
				for (int ed = 0; ed < n_ed; ++ed)
				{
					settings.EdgeDetection = edge_detects[ed];
					weights.SetSettings(settings);
					weights.CalculateRegion(cx, cy, min_room);
					_snprintf(name, sizeof(name) / sizeof(name[0]), "%s-%s-%s-%s", px_red_strs[pr], noise_red_strs[nr], accentuate_strs[a], edge_detect_strs[ed]);
					weights.SaveImage(name);
				}
			}
		}
	}
}
#endif
