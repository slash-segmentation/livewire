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

const QBrush QtPanel::NotCalculatedBrush(QColor(0xFF, 0, 0, 0x30));
const QBrush QtPanel::CalculatingBrush(QColor(0xFF, 0xFF, 0, 0x30));

QtPanel::QtPanel(QWidget *parent) : QWidget(parent),
	livewire(NULL), wrapwire(NULL), w(0), h(0), showCalculationBlocks(false), availSize(1024)
{
	this->setMouseTracking(true);
	this->setWindowTitle("Livewire Demo");
	this->resize(200, 150);
	this->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(this, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(ShowContextMenu(const QPoint&)));
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
}

void QtPanel::ShowContextMenu(const QPoint &p)
{
	QPoint globalPos = this->mapToGlobal(p);

	QMenu menu(this);

	QAction setImage("Set &Image...", this);

	QAction showCalcBoxes("Show Calculation &Boxes", this); showCalcBoxes.setCheckable(true); showCalcBoxes.setChecked(this->showCalculationBlocks);

#define CHECK_BOX_OPT(name, text, checked, menu) QAction name(text, this); name.setCheckable(true); name.setChecked(checked); menu.addAction(&name);

	QMenu availSizes("&Minimum Window", this);
	CHECK_BOX_OPT(as512 , "&512x512",   this->availSize ==  512, availSizes);
	CHECK_BOX_OPT(as1024, "&1024x1024", this->availSize == 1024, availSizes);
	CHECK_BOX_OPT(as2048, "&2048x2048", this->availSize == 2048, availSizes);
	CHECK_BOX_OPT(as4096, "&4096x4096", this->availSize == 4096, availSizes);

	Weights::Settings settings = this->livewire->GetSettings();

	QAction grysclSettings("Default &Grayscale Settings", this); grysclSettings.setCheckable(true); grysclSettings.setChecked(settings == Weights::GrayscaleSettings);
	QAction colorSettings ("Default &Color Settings",     this); colorSettings.setCheckable(true);  colorSettings.setChecked (settings == Weights::ColorSettings    );

	QMenu sttngs("&Settings", this);
	QMenu method("&Coalescing Method", this); sttngs.addMenu(&method);
	CHECK_BOX_OPT(red,      "&Red Channel",     settings.Method == Weights::RedChannel,   method);
	CHECK_BOX_OPT(green,    "&Green Channel",   settings.Method == Weights::GreenChannel, method);
	CHECK_BOX_OPT(blue,     "&Blue Channel",    settings.Method == Weights::BlueChannel,  method);
	CHECK_BOX_OPT(avg,      "&Average RGB",     settings.Method == Weights::AvgRGB,       method);
	CHECK_BOX_OPT(luma,     "&Luma",            settings.Method == Weights::Luma,         method);
	CHECK_BOX_OPT(luma601,  "Luma Rec &601",    settings.Method == Weights::Luma601,      method);
	CHECK_BOX_OPT(lumaSMPTE,"Luma SMPTE &240M", settings.Method == Weights::LumaSMPTE,    method);
	CHECK_BOX_OPT(hsv,      "Weighted HS&V",    settings.Method == Weights::WeightedHSV,  method);
	CHECK_BOX_OPT(hsl,      "Weighted &HSL",    settings.Method == Weights::WeightedHSL,  method);
	CHECK_BOX_OPT(hsi,      "Weighted HS&I",    settings.Method == Weights::WeightedHSI,  method);

	QMenu pixel("&Pixel Reduction", this); sttngs.addMenu(&pixel);
	CHECK_BOX_OPT(pixelNone,    "&None",                  settings.PixelReduction == Weights::NoPixelReduction,  pixel);
	CHECK_BOX_OPT(pixelMedian2, "&Median (2px window)",   settings.PixelReduction == Weights::Median2pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian3, "&Median (3px window)",   settings.PixelReduction == Weights::Median3pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian4, "&Median (4px window)",   settings.PixelReduction == Weights::Median4pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian5, "&Median (5px window)",   settings.PixelReduction == Weights::Median5pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMean2,   "&Mean (2px window)",     settings.PixelReduction == Weights::Mean2pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean3,   "&Mean (3px window)",     settings.PixelReduction == Weights::Mean3pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean4,   "&Mean (4px window)",     settings.PixelReduction == Weights::Mean4pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean5,   "&Mean (5px window)",     settings.PixelReduction == Weights::Mean5pxWindow,     pixel);
	CHECK_BOX_OPT(pixelGaus3,   "&Gaussian (3px window)", settings.PixelReduction == Weights::Gaussian3pxWindow, pixel);
	CHECK_BOX_OPT(pixelGaus4,   "&Gaussian (4px window)", settings.PixelReduction == Weights::Gaussian4pxWindow, pixel);
	CHECK_BOX_OPT(pixelGaus5,   "&Gaussian (5px window)", settings.PixelReduction == Weights::Gaussian5pxWindow, pixel);

	QMenu noise("&Noise Reduction", this); sttngs.addMenu(&noise);
	CHECK_BOX_OPT(noiseNone,    "&None",                  settings.NoiseReduction == Weights::NoNoiseReduction,        noise);
	CHECK_BOX_OPT(noiseMedian3, "&Median (3px window)",   settings.NoiseReduction == Weights::MedianFilter3pxWindow,   noise);
	CHECK_BOX_OPT(noiseMedian5, "&Median (5px window)",   settings.NoiseReduction == Weights::MedianFilter5pxWindow,   noise);
	CHECK_BOX_OPT(noiseMean3,   "&Mean (3px window)",     settings.NoiseReduction == Weights::MeanFilter3pxWindow,     noise);
	CHECK_BOX_OPT(noiseMean5,   "&Mean (5px window)",     settings.NoiseReduction == Weights::MeanFilter5pxWindow,     noise);
	CHECK_BOX_OPT(noiseGaus3,   "&Gaussian (3px window)", settings.NoiseReduction == Weights::GaussianFilter3pxWindow, noise);
	CHECK_BOX_OPT(noiseGaus5,   "&Gaussian (5px window)", settings.NoiseReduction == Weights::GaussianFilter5pxWindow, noise);

	QMenu edged("&Edge Detection", this); sttngs.addMenu(&edged);
	CHECK_BOX_OPT(edgeNone,  "&None",  settings.EdgeDetection == Weights::NoEdgeDetection, edged);
	CHECK_BOX_OPT(edgeSobel, "&Sobel", settings.EdgeDetection == Weights::Sobel,           edged);

	QMenu accnt("&Accentuation", this); sttngs.addMenu(&accnt);
	CHECK_BOX_OPT(accntNone,    "&None",    settings.Accentuation == Weights::NoAccentuation, accnt);
	CHECK_BOX_OPT(accntSigmoid, "&Sigmoid", settings.Accentuation == Weights::Sigmoid,        accnt);

	QAction invert("&Invert", this); invert.setCheckable(true); invert.setChecked(settings.Invert); sttngs.addAction(&invert);

	menu.addAction(&setImage);
	menu.addAction(&showCalcBoxes);
	menu.addMenu(&availSizes);
	menu.addSeparator();
	menu.addAction(&grysclSettings);
	menu.addAction(&colorSettings);
	menu.addMenu(&sttngs);

	QAction* selectedItem = menu.exec(globalPos);
	if (selectedItem)
	{
		if (selectedItem == &setImage)
		{
			QString fileName = QFileDialog::getOpenFileName(this, "Open Image", QDir::currentPath(), "Image Files (*.png *.jpg *.bmp)");
			if (!fileName.isNull()) { QImage img(fileName); this->SetImage(img); }
		}
		else if (selectedItem == &showCalcBoxes) { this->showCalculationBlocks = !this->showCalculationBlocks; this->update(); }
		else if (selectedItem == &as512 ) { this->availSize =  512; }
		else if (selectedItem == &as1024) { this->availSize = 1024; }
		else if (selectedItem == &as2048) { this->availSize = 2048; }
		else if (selectedItem == &as4096) { this->availSize = 4096; }
/*		else if (selectedItem == &grysclSettings) { this->weights->ChangeSettings(Weights::GrayscaleSettings); }
		else if (selectedItem == &colorSettings)  { this->weights->ChangeSettings(Weights::ColorSettings    ); }
		else if (selectedItem == &red)       { settings.Method = Weights::RedChannel;   this->weights->ChangeSettings(settings); }
		else if (selectedItem == &green)     { settings.Method = Weights::GreenChannel; this->weights->ChangeSettings(settings); }
		else if (selectedItem == &blue)      { settings.Method = Weights::BlueChannel;  this->weights->ChangeSettings(settings); }
		else if (selectedItem == &avg)       { settings.Method = Weights::AvgRGB;       this->weights->ChangeSettings(settings); }
		else if (selectedItem == &luma)      { settings.Method = Weights::Luma;         this->weights->ChangeSettings(settings); }
		else if (selectedItem == &luma601)   { settings.Method = Weights::Luma601;      this->weights->ChangeSettings(settings); }
		else if (selectedItem == &lumaSMPTE) { settings.Method = Weights::LumaSMPTE;    this->weights->ChangeSettings(settings); }
		else if (selectedItem == &hsv)       { settings.Method = Weights::WeightedHSV;  this->weights->ChangeSettings(settings); }
		else if (selectedItem == &hsl)       { settings.Method = Weights::WeightedHSL;  this->weights->ChangeSettings(settings); }
		else if (selectedItem == &hsi)       { settings.Method = Weights::WeightedHSI;  this->weights->ChangeSettings(settings); }*/

/*	CHECK_BOX_OPT(pixelNone,    "&None",                  settings.PixelReduction == Weights::NoPixelReduction,  pixel);
	CHECK_BOX_OPT(pixelMedian2, "&Median (2px window)",   settings.PixelReduction == Weights::Median2pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian3, "&Median (3px window)",   settings.PixelReduction == Weights::Median3pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian4, "&Median (4px window)",   settings.PixelReduction == Weights::Median4pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMedian5, "&Median (5px window)",   settings.PixelReduction == Weights::Median5pxWindow,   pixel);
	CHECK_BOX_OPT(pixelMean2,   "&Mean (2px window)",     settings.PixelReduction == Weights::Mean2pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean3,   "&Mean (3px window)",     settings.PixelReduction == Weights::Mean3pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean4,   "&Mean (4px window)",     settings.PixelReduction == Weights::Mean4pxWindow,     pixel);
	CHECK_BOX_OPT(pixelMean5,   "&Mean (5px window)",     settings.PixelReduction == Weights::Mean5pxWindow,     pixel);
	CHECK_BOX_OPT(pixelGaus3,   "&Gaussian (3px window)", settings.PixelReduction == Weights::Gaussian3pxWindow, pixel);
	CHECK_BOX_OPT(pixelGaus4,   "&Gaussian (4px window)", settings.PixelReduction == Weights::Gaussian4pxWindow, pixel);
	CHECK_BOX_OPT(pixelGaus5,   "&Gaussian (5px window)", settings.PixelReduction == Weights::Gaussian5pxWindow, pixel);

	CHECK_BOX_OPT(noiseNone,    "&None",                  settings.NoiseReduction == Weights::NoNoiseReduction,        noise);
	CHECK_BOX_OPT(noiseMedian3, "&Median (3px window)",   settings.NoiseReduction == Weights::MedianFilter3pxWindow,   noise);
	CHECK_BOX_OPT(noiseMedian5, "&Median (5px window)",   settings.NoiseReduction == Weights::MedianFilter5pxWindow,   noise);
	CHECK_BOX_OPT(noiseMean3,   "&Mean (3px window)",     settings.NoiseReduction == Weights::MeanFilter3pxWindow,     noise);
	CHECK_BOX_OPT(noiseMean5,   "&Mean (5px window)",     settings.NoiseReduction == Weights::MeanFilter5pxWindow,     noise);
	CHECK_BOX_OPT(noiseGaus3,   "&Gaussian (3px window)", settings.NoiseReduction == Weights::GaussianFilter3pxWindow, noise);
	CHECK_BOX_OPT(noiseGaus5,   "&Gaussian (5px window)", settings.NoiseReduction == Weights::GaussianFilter5pxWindow, noise);
	*/
/*		else if (selectedItem == &edgeNone)     { settings.EdgeDetection = Weights::NoEdgeDetection; this->weights->ChangeSettings(settings); }
		else if (selectedItem == &edgeSobel)    { settings.EdgeDetection = Weights::Sobel;           this->weights->ChangeSettings(settings); }
		else if (selectedItem == &accntNone)    { settings.Accentuation = Weights::NoAccentuation;   this->weights->ChangeSettings(settings); }
		else if (selectedItem == &accntSigmoid) { settings.Accentuation = Weights::Sigmoid;          this->weights->ChangeSettings(settings); }
		else if (selectedItem == &invert)       { settings.Invert = !settings.Invert;                this->weights->ChangeSettings(settings); }*/
		else { /* unknown item was chosen */ }
	}
	else { /* nothing was chosen */ }
}


void QtPanel::SetImage(QImage &img)
{
	this->Cleanup();
	if (!img.isNull() && img.width() > 1 && img.height() > 1)
	{
		// Basic image updating
		bool gray = img.format() == QImage::Format_Indexed8 && img.isGrayscale();
		this->image = gray ? img.copy() : img.convertToFormat(QImage::Format_RGB32);
		//this->image.SetResolution(this->dpiX, this->dpiY); // TODO: DPI stuff
		uint w = this->w = img.width(), h = this->h = img.height();

		// Setup the calculators
		AddEventsToThreaded(this->livewire = new LivewireCalculator());
		this->livewire->SetSettings(img.allGray() ? Weights::GrayscaleSettings : Weights::ColorSettings);
		this->livewire->SetImage(this->image.bits(), w, h, gray ? Weights::GrayscaleByte : Weights::RGB, this->image.bytesPerLine());

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
		if (this->showCalculationBlocks)
		{
			/*TODO: g.setPen(Qt::black);
			QVector<QPoint> done, doing, xdone;
			uint block_size = this->weights->GetBlocksCalculated(done, doing, xdone);
			g.setBrush(NotCalculatedBrush);
			for (int i = 0; i < xdone.count(); ++i)	g.drawRect(xdone[i].x(), xdone[i].y(), block_size, block_size);
			g.setBrush(CalculatingBrush);
			for (int i = 0; i < doing.count(); ++i)	g.drawRect(doing[i].x(), doing[i].y(), block_size, block_size);*/
		}

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
			// TODO: Start wrapping from a new point
			//this->wrapwire = lw;
			//AddEventsToThreaded(this->livewire = new LivewireCalculator(this->weights));
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
		this->livewire->Start(p.x(), p.y(), this->availSize / 2);

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
