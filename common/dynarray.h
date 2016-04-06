#ifndef RE_DYNARRAY_H
#define RE_DYNARRAY_H

#include "assert.h"
#include <cstring> // memcpy
#include <cstdlib> // malloc
#include <stddef.h> // ptrdiff

namespace Re {

// Equivalent class as std::vector (calling it "vector" was a great mistake...)

template<typename T>
class DynArray
{
private:
    static const int MinimalAllocSize = 8;


public:

    typedef T value_type;
    typedef int size_type;

    value_type* mData = 0 ;
    size_type mSize = 0;
    size_type mCapacity = 0;

    typedef value_type*         iterator;
    typedef const value_type*   const_iterator;

    ~DynArray() {

        if (mData) {
            free(mData);
            mData = 0;
        }
    }

    inline iterator begin() {
        return mData;
    }

    inline iterator end() {
        return mData + mSize;
    }

    inline const_iterator begin() const {
        return mData;
    }

    inline const_iterator end() const {
        return mData + mSize;
    }

    inline bool empty() const {
        return mSize == 0;
    }

    inline size_type size() const{
        return mSize;
    }

    inline size_type capacity() const {
        return mCapacity;
    }

    inline value_type& operator[](int i)
    {
        assert((i >= 0 && i < mSize));

        return mData[i];
    }

    inline const value_type& operator[](int i) const {
        assert((i >= 0 && i < mSize));

        return mData[i];
    }

    inline void clear() {
        if (mData) {
            mSize = mCapacity = 0;
            free(mData);
            mData = NULL;
        }
    }

    inline value_type& front() {
        assert(mSize > 0);

        return mData[0];
    }

    inline const value_type& front() const {
        assert(mSize > 0);

        return mData[0];
    }

    inline value_type& back() {
        assert(mSize > 0);

        return mData[mSize - 1];
    }

    inline const value_type& back() const {
        assert(mSize > 0);

        return mData[mSize - 1];
    }

    inline void swap(DynArray& rhs) {

        if (&rhs != this) {

            int rhsSize = rhs.mSize;
            int rhsCapacity = rhs.mCapacity;
            value_type* rhsData = rhs.mData;

            rhs.mSize = mSize;
            mSize = rhsSize;

            rhs.mCapacity = mCapacity;
            mCapacity = rhsCapacity;

            rhs.Data = mData;
            mData = rhsData;
        }
    }


    // affects the size
    inline void resize(size_type neededSize) {

        if (neededSize == mSize) {
            return; // Do nothing
        }

        if (neededSize > mCapacity) {
            size_type neededCapacity = needed_capacity(neededSize);
            reserve(neededCapacity);
        } else {
            mSize = neededSize;
        }
    }

    //
    // Increase the capacity of the container to a value that's greater or equal to new_cap.
    // If new_cap is greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
    inline void reserve(size_type new_cap)
    {
        if (new_cap <= mCapacity) {
            return; // Do nothing
        }
        value_type* newData = (value_type*)malloc((size_t)new_cap * sizeof(value_type));
        memcpy(newData, mData, (size_t)mSize * sizeof(value_type));
        free(mData);
        mData = newData;
        mCapacity = new_cap;
    }

    inline void push_back(const value_type& value) {
        if (mSize == mCapacity) {
            size_type neededSize = mSize + 1;
            reserve(needed_capacity(neededSize));
        }
        mData[mSize++] = value;
    }

    inline void pop_back() {
        assert(mSize > 0);
        --mSize;
    }

    inline const_iterator erase(const_iterator it) {
        assert(it >= mData && it < mData + mSize);



            const ptrdiff_t off = it - mData;
            memmove(mData + off, mData + off + 1, ((size_t)mSize - (size_t)off - 1) * sizeof(value_type));
            mSize--;
            return mData + off;

            return erase(it, it + 1);
    }

    const_iterator erase(const_iterator first, const_iterator last) {

        const size_type first_index = size_type(first - mData);
        const size_type last_index = size_type(last - mData);
        const size_type size_removed = last_index - first_index;
        const size_type size_to_move = mSize - last_index;

        memmove(mData + first_index, mData + last_index, size_to_move * sizeof(value_type));

        mSize -= size_removed;
        return mData + first_index;
    }

    //inline iterator             insert(const_iterator it, const value_type& v)  { IM_ASSERT(it >= Data && it <= Data+Size); const ptrdiff_t off = it - Data; if (Size == Capacity) reserve(Capacity ? Capacity * 2 : 4); if (off < (int)Size) memmove(Data + off + 1, Data + off, ((size_t)Size - (size_t)off) * sizeof(value_type)); Data[off] = v; Size++; return Data + off; }

private:
    inline int needed_capacity(int neededSize) {

        // If the array is not initialize, set the minimal size.
        // Otherwise increase the size by 50%
        int newCapacity = !mCapacity ? MinimalAllocSize : (mCapacity + (mCapacity / 2));

        // Use the greatest of both neededSize and newCapacity
        return newCapacity > neededSize ? newCapacity : neededSize;
    }
};

} // namespace Re


#endif // RE_DYNARRAY_H
