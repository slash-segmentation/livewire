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

#ifndef LIVEWIREQTPANEL_H
#define LIVEWIREQTPANEL_H

#include "LivewireCalculator.h"

#include <QWidget>
#include <QImage>
#include <QPoint>
#include <QVector>

namespace Livewire
{
	class QtPanel : public QWidget
	{
		Q_OBJECT

	private:
		/// <summary>The pen used to draw the active livewire</summary>
		static const QPen LivewirePen;

		/// <summary>The pen used to draw the static previous wires</summary>
		static const QPen PrevwirePen;

		/// <summary>The pen used to draw the wrapping livewire</summary>
		static const QPen WrapwirePen;

		static const QBrush MouseBrush;
		static const QBrush PointBrush;

		static const QBrush NotCalculatedBrush;
		static const QBrush CalculatingBrush;

	private:
		QImage image;
		uint w, h;
		QPoint mouse;
		QVector<QPoint> points;
		LivewireCalculator *livewire;
		LivewireCalculator *wrapwire;
		QImage prevwires;

		bool showCalculationBlocks;
		uint availSize;

		void Cleanup();

	public:
		QtPanel(QWidget *parent = NULL);
		~QtPanel();

		void SetImage(QImage &img);
		/// <summary>Reset the current state of the livewires</summary>
		void Reset();

	protected:
		void paintEvent(QPaintEvent *evnt);
		void mouseMoveEvent(QMouseEvent *evnt);
		void mousePressEvent(QMouseEvent *evnt);
		void mouseDoubleClickEvent(QMouseEvent *evnt);

	private:
		/// <summary>Handles mouse clicking and double clicking. </summary>
		/// <param name="evnt">The mouse event arguments</param>
		/// <param name="dbl">True if the mouse double clicked</param>
		void MouseClicked(QMouseEvent *evnt, bool dbl);

		void AddEventsToThreaded(const Threaded *t) const;

		/// <summary>Moves the point to be within the bounds of the image</summary>
		/// <param name="p">The point to move</param>
		/// <returns>The adjusted point</returns>
		void getPointInBounds(QPoint &p);

	private slots:
		void ShowContextMenu(const QPoint &p);
	};
}

#endif
