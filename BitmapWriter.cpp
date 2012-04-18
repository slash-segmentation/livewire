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

void Livewire::WriteBitmap(const uint *data, uint w, uint h, const char *name)
{
	QImage b(w, h, QImage::Format_Indexed8);
	b.setColorTable(*GetGrayscaleColorTable());
	uint max = 0, wh = w*h;
	for (uint i = 0; i < wh; ++i)
		if (data[i] > max)
			max = data[i];
	for (uint y = 0, off = 0; y < h; ++y, off += w)
	{
		byte *line = b.scanLine(y);
		for (uint x = off, X = 0; X < w; ++x, ++X)
			line[X] = (byte)(data[x] * 255 / max);
	}
	b.save(QString::number(QDateTime::currentMSecsSinceEpoch()) + "-" + QString(name) + ".png");
}

void Livewire::WriteBitmap(const byte *data, uint w, uint h, const char *name)
{
	QImage b(w, h, QImage::Format_Indexed8);
	b.setColorTable(*GetGrayscaleColorTable());
	for (uint y = 0, off = 0; y < h; ++y, off += w)
		memcpy(b.scanLine(y), data+off, w);
	b.save(QString::number(QDateTime::currentMSecsSinceEpoch()) + "-" + QString(name) + ".png");
}
