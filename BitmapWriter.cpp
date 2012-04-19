/**

Livewire - Core code for running the Livewire algorithm
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

// This class is optional, and only used if SAVE_WEIGHT_IMAGE or SAVE_SCORE_IMAGE are defined

#include "BitmapWriter.h"

#include <QImage>
#include <QVector>
#include <QDateTime>

static QVector<QRgb> *grayscale_color_table;

static const QVector<QRgb> *GetGrayscaleColorTable()
{
	if (grayscale_color_table == NULL)
	{
		grayscale_color_table = new QVector<QRgb>();
		grayscale_color_table->reserve(256);
		for (int i = 0; i < 256; ++i)
			grayscale_color_table->append(qRgb(i, i, i));
		grayscale_color_table->squeeze();
	}
	return grayscale_color_table;
}

void Livewire::WriteBitmap(const SparseMatrix<uint> &data, uint w, uint h, const char *name)
{
	QImage b(w, h, QImage::Format_Indexed8);
	b.setColorTable(*GetGrayscaleColorTable());
	uint max = 0;
	for (uint y = 0; y < h; ++y)
	{
		for (uint x = 0; x < w; ++x)
		{
			uint val = data.Get(x, y);
			if (val > max)
				max = val;
		}
	}
	for (uint y = 0; y < h; ++y)
	{
		byte *line = b.scanLine(y);
		for (uint x = 0; x < w; ++x)
			line[x] = (byte)(data.Get(x, y) * 255 / max);
	}
	b.save(QString::number(QDateTime::currentMSecsSinceEpoch()) + "-" + QString(name) + ".png");
}

void Livewire::WriteBitmap(const RawSparseMatrix data, uint w, uint h, uint bs, const char *name)
{
	QImage b(w, h, QImage::Format_Indexed8);
	b.setColorTable(*GetGrayscaleColorTable());
	const uint W = ScaleBack(w, bs), H = ScaleBack(h, bs), W1 = W - 1, H1 = H - 1, bw_ = w % bs, bh_ = h % bs;

	for (uint Y = 0, y = 0, I = 0; Y < H; ++Y, y += bs)
	{
		const uint bh = (Y == H1) ? bh_ : bs;
		for (uint X = 0, x = 0; X < W; ++X, x += bs, ++I)
		{
			const uint bw = (X == W1) ? bw_ : bs;
			if (data[I])	for (uint y_ = 0; y_ < bh; ++y_) { memcpy(b.scanLine(y_ + y) + x, data[I]+y_*bw, bw); }
			else			for (uint y_ = 0; y_ < bh; ++y_) { memset(b.scanLine(y_ + y) + x, 0,             bw); }
		}
	}
	b.save(QString::number(QDateTime::currentMSecsSinceEpoch()) + "-" + QString(name) + ".png");
}
