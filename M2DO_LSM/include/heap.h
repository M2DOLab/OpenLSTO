/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+lsm@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// Adapted from Scikit-FMM: https://github.com/scikit-fmm/scikit-fmm

#ifndef _HEAP_H
#define _HEAP_H

/*! \file Heap.h
    \brief An implementation of a heap data structure (binary tree).
 */

/*! \brief An implementation of a heap data structure (binary tree).

    The heap is stored as a contiguous vector array. The children of a given
    entry i are indexed as follows:

        2*i + 1 = Left-hand child
        2*i + 2 = Right-hand child

    The heap is used to maintain an ascending order priority queue of unsigned
    node distances from the zero iso-contour of the level set.

    For efficiency, the heap is stored as a contiguous std::vector array,
    rather than a linked list. The constructor needs to know the maximum number
    of entries that will be added to the heap. This will number be computed in
    the FastMarchingMethod::initialiseHeap method (which calculates the number
    of far field nodes).
 */
class Heap
{
public:
    //! Constructor.
    /*! \param maxLength_
            The maximum number of entries in the heap.

        \param isTest_
            Whether to perform self testing of the heap (optional).
     */
    Heap(unsigned int, bool isTest_ = false);

    //! Push a value onto the heap.
    /*! \param address_
            The address of the element to push (its array index).

        \param value
            The value of the element.

        \return
            The index of the value in the heap.
     */
    unsigned int push(unsigned int, double);

    //! Pop the top value from the heap.
    /*! \param address_
            The address (index) of top entry in the heap.

        \param value
            The value of the top entry in the heap.
     */
    void pop(unsigned int&, double&);

    //! Set a specific heap entry.
    /*! \param index
            The index in the heap.

        \param newDistance
            The new distance value to insert into the heap.
     */
    void set(unsigned int, double);

    //! Test whether the heap is empty.
    /*! \return
            Whether the heap is empty (true) or contains entries (false).
     */
    bool empty() const;

    //! Return the value of the top entry in the heap.
    /*! \return
            The value currently at the top of the heap.
     */
    const double& peek() const;

    //! Return the current size of the heap.
    /* \return
            The current size of the heap.
     */
    const unsigned int& size() const;

private:
    //! Test that the heap is correct.
    void test() const;

    //! Sift a value up the heap.
    /*! \param pos
            The position in the heap.
     */
    void siftUp(unsigned int);

    //! Sift a value down the heap.
    /*! \param startPos
            The starting position in the heap.

        \param pos
            The current position in the heap.
     */
    void siftDown(unsigned int, unsigned int);

    /// The maximum number of entries in the heap.
    unsigned int maxLength;

    /// The current size of the heap.
    unsigned int heapLength;

    /// The current size of the list.
    unsigned int listLength;

    /// The unsigned distance from the zero level set iso-contour.
    std::vector<double> distance;

    /// The index of entries into the distance (or address) array.
    std::vector<unsigned int> heap;

    /// The (original) grid address of each element in the heap.
    std::vector<unsigned int> address;

    /// A map from the index of distance (or address) to the current location of the heap element.
    std::vector<unsigned int> backPointer;

    /// Whether self testing of the heap is active.
    bool isTest;
};

#include "../src/heap.cpp"

#endif  /* _HEAP_H */
