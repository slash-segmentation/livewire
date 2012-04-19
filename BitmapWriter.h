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

#ifndef BITMAP_WRITER_H
#define BITMAP_WRITER_H

#include "general.h"

namespace Livewire
{
	/// <summary>Debug utility to save a sparse matrix as a grayscale image</summary>
	/// <param name="data">The data to save</param>
	/// <param name="w">The width of the data</param>
	/// <param name="h">The height of the data</param>
	/// <param name="name">The name of the data</param>
	void WriteBitmap(const SparseMatrix<uint> &data, uint w, uint h, const char *name);

	typedef byte ** RawSparseMatrix;

	/// <summary>Debug utility to save a raw sparse matrix as a grayscale image</summary>
	/// <param name="data">The data to save</param>
	/// <param name="w">The width of the data</param>
	/// <param name="h">The height of the data</param>
	/// <param name="bs">The block-size of the data</param>
	/// <param name="name">The name of the data</param>
	void WriteBitmap(const RawSparseMatrix data, uint w, uint h, uint bs, const char *name);
}

#endif
