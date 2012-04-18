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

#include <QtGui>

using namespace Livewire;

const QPen QtPanel::LivewirePen(Qt::red, 1.5f);
const QPen QtPanel::PrevwirePen(Qt::yellow, 1.5f);
const QPen QtPanel::WrapwirePen(QColor(255, 128, 128), 1.5f);

const QBrush QtPanel::MouseBrush(Qt::red);
const QBrush QtPanel::PointBrush(Qt::yellow);

QtPanel::QtPanel(QWidget *parent) : QWidget(parent),
	weights(NULL), livewire(NULL), wrapwire(NULL), w(0), h(0)
{
	this->setMouseTracking(true);
	this->setWindowTitle("Livewire Demo");
	this->resize(200, 150);
}
QtPanel::~QtPanel() { this->Cleanup(); }
void QtPanel::Cleanup()
{
	if (this->livewire != NULL)
	{
		this->livewire->Stop();
		delete this->livewire;
		this->livewire = NULL;
		this->prevwires = QImage();
	}
	if (this->wrapwire != NULL)
	{
		this->wrapwire->Stop();
		delete this->wrapwire;
		this->wrapwire = NULL;
	}
	if (this->weights != NULL)
	{
		this->weights->Stop();
		delete this->weights;
		this->weights = NULL;
	}
}

void QtPanel::SetImage(QImage &img)
{
	this->Cleanup();
	if (!img.isNull() && img.width() > 1 && img.height() > 1)
	{
		// Basic image updating
		bool gray = img.format() == QImage::Format_Indexed8 && img.isGrayscale();
		this->image = gray ? img : img.convertToFormat(QImage::Format_RGB32);
		//this->image.SetResolution(this->dpiX, this->dpiY); // TODO: DPI stuff
		uint w = this->w = img.width(), h = this->h = img.height();

		// Setup the calculators
		this->weights = new WeightCalculator(w, h, WeightCalculator::GrayScaleSettings);
		this->weights->SetImage(this->image.bits(), gray ? WeightCalculator::GrayscaleByte : WeightCalculator::RGB, this->image.bytesPerLine());
		AddEventsToThreaded(this->livewire = new LivewireCalculator(this->weights));

		// Update variables for the new image size
		this->prevwires = QImage(w, h, QImage::Format_ARGB32_Premultiplied);
		this->prevwires.fill(Qt::transparent);
		this->resize(this->image.size());
	}
	else
	{
		// Clear variables
		this->image = QImage();
		this->w = 0;
		this->h = 0;
		this->resize(200, 150);
	}

	// Update general variables
	this->points.clear();

	this->update();
}

void QtPanel::paintEvent(QPaintEvent *evnt)
{
	(evnt); // unreferenced

	QPainter g(this);
	g.setRenderHint(QPainter::Antialiasing);
	// TODO: dpi stuff?
	//g.Clear(this->BackColor);

	if (!this->image.isNull())
	{
		// The weights are calculated so we can do normal painting
		g.drawImage(0, 0, this->image); // draw the image itself

		// Draw the not calculated / calculating blocks
		g.setPen(Qt::black);
		QVector<QPoint> done, doing, not_done;
		uint block_size = this->weights->GetBlocksCalculated(done, doing, not_done);
		g.setBrush(QBrush(QColor(0xFF, 0, 0, 0x30)));
		for (int i = 0; i < not_done.count(); ++i)
			g.drawRect(not_done[i].x(), not_done[i].y(), block_size, block_size);
		g.setBrush(QBrush(QColor(0xFF, 0xFF, 0, 0x30)));
		for (int i = 0; i < doing.count(); ++i)
			g.drawRect(doing[i].x(), doing[i].y(), block_size, block_size);

		if (this->livewire != NULL)
		{
			// Draw the previous, wrapping, and current livewires
			g.drawImage(0, 0, this->prevwires);
			if (this->points.size() > 0)
			{
				g.setBrush(Qt::NoBrush);
				if (this->wrapwire != NULL)
				{
					g.setPen(WrapwirePen);
					this->wrapwire->DrawTrace(this->mouse.x(), this->mouse.y(), g);
				}
				g.setPen(LivewirePen);
				this->livewire->DrawTrace(this->mouse.x(), this->mouse.y(), g);
			}
		}

		g.setPen(Qt::NoPen); g.setBrush(PointBrush);
		for (int i = 0; i < this->points.size(); ++i)
			g.drawEllipse(this->points[i], 4, 4); // draw the seed points

		if (!this->mouse.isNull()) { g.setBrush(MouseBrush); g.drawEllipse(this->mouse, 4, 4); } // draw the mouse indicator

		if (this->livewire->IsExecuting())
		{
			// The livewire data is not completely calculated, so we draw a progress message
			float percent = (float)this->livewire->GetProgress() / this->livewire->GetTotalProgress();
			QString text = "Calculating livewire... " + QString::number((int)(percent * 100)) + "%";
			QFontMetrics fm = g.fontMetrics();
			QSize sz = fm.size(Qt::TextSingleLine, text);
			g.setPen(Qt::black); g.setBrush(Qt::white);   g.drawRect(0, 0, sz.width()+4, sz.height()+4);
			g.setPen(Qt::NoPen); g.setBrush(Qt::blue);    g.drawRect(1, 1, (sz.width()+2) * percent, sz.height()+2);
			g.setPen(Qt::black); g.setBrush(Qt::NoBrush); g.drawText(2, 2+fm.ascent(), text);
		}
	}
}

void QtPanel::mouseMoveEvent(QMouseEvent *evnt)
{
	QPoint p = evnt->pos();
	this->getPointInBounds(p);
	if (this->mouse != p)
	{
		this->mouse = p;
		this->update();
	}
}
void QtPanel::mousePressEvent(QMouseEvent *evnt)       { this->MouseClicked(evnt, false); }
void QtPanel::mouseDoubleClickEvent(QMouseEvent *evnt) { this->MouseClicked(evnt, true); }
void QtPanel::MouseClicked(QMouseEvent *evnt, bool dbl)
{
	if (evnt->button() == Qt::LeftButton && this->livewire != NULL)
	{
		bool have_points = this->points.size() > 0;
		QPoint p = evnt->pos();
		this->getPointInBounds(p);
		this->mouse = p;

		LivewireCalculator *lw = this->livewire;

		if (dbl && this->wrapwire != NULL)
		{
			// On double click draw the wrapping livewire
			QPainter g(&this->prevwires);
			g.setRenderHint(QPainter::Antialiasing);
			g.setPen(PrevwirePen);
			g.setBrush(Qt::NoBrush);
			this->wrapwire->DrawTrace(p.x(), p.y(), g);
			this->wrapwire = NULL;
		}
		else if (have_points && this->wrapwire == NULL)
		{
			// Start wrapping from a new point
			this->wrapwire = lw;
			AddEventsToThreaded(this->livewire = new LivewireCalculator(this->weights));
		}

		if (have_points)
		{
			// On a single click save the livewire
			{
				QPainter g(&this->prevwires);
				g.setRenderHint(QPainter::Antialiasing);
				g.setPen(PrevwirePen);
				g.setBrush(Qt::NoBrush);
				lw->DrawTrace(p.x(), p.y(), g);
			}
		}

		// Start preparing the next livewire
		this->livewire->Start(p.x(), p.y(), 512);

		this->points.push_back(p);
		this->update();
	}
}


void QtPanel::AddEventsToThreaded(const Threaded *t) const
{
	QObject::connect(t, SIGNAL(ProgressChanged(int)), this, SLOT(update()));
	QObject::connect(t, SIGNAL(finished()),           this, SLOT(update()));
}

void QtPanel::getPointInBounds(QPoint &p)
{
	if (p.x() < 0) { p.setX(0); } else if (p.x() >= (int)this->w) { p.setX(this->w - 1); }
	if (p.y() < 0) { p.setY(0); } else if (p.y() >= (int)this->h) { p.setY(this->h - 1); }
}

void QtPanel::Reset()
{
	this->points.clear();
	if (this->livewire != NULL)
	{
		this->livewire->Stop();
		this->prevwires.fill(Qt::transparent);
	}
	if (this->wrapwire != NULL)
	{
		this->wrapwire->Stop();
		this->wrapwire = NULL;
	}
	this->update();
}
