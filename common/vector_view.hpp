#pragma once

#include <vector>

namespace rt {
	template <class T>
	class vector_view {
	public:
		vector_view(T *ptr, std::size_t size) :_ptr(ptr), _size(size) {}

	private:
		T *_ptr;
		std::size_t _size;
	};
}