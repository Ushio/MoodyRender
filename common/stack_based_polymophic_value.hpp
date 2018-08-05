#pragma once

#include <type_traits>
#include <algorithm>
#include <string>

namespace rt {
	namespace sbpv_details {
		template <class...>
		struct is_one_of {
			static constexpr bool value = false;
		};
		template <class S, class T, class... Ts>
		struct is_one_of<S, T, Ts...> {
			static constexpr bool value = std::is_same<S, T>::value || is_one_of<S, Ts...>::value;
		};

		template <class...>
		struct is_base_all_of {
			static constexpr bool value = true;
		};
		template <class TBase, class T, class... Ts>
		struct is_base_all_of<TBase, T, Ts...> {
			static constexpr bool value = std::is_base_of<TBase, T>::value && is_base_all_of<TBase, Ts...>::value;
		};

		template <class TBase>
		class InstanceManager {
		public:
			virtual TBase *ptr(void *p) const = 0;
			virtual const TBase *ptr(const void *p) const = 0;
			virtual TBase *copyConstruct(void *p, const void *src) const = 0;
			virtual void destruct(void *p) const = 0;
		};

		template <class TBase, class T>
		class InstanceManagerT : public InstanceManager<TBase> {
		public:
			static_assert(std::is_base_of<TBase, T>::value, "T is must be inherits TBase");

			TBase *ptr(void *p) const override {
				return static_cast<T *>(p);
			}
			virtual const TBase *ptr(const void *p) const {
				return static_cast<const T *>(p);
			}
			TBase *copyConstruct(void *p, const void *src) const override {
				const T *srcDerived = static_cast<const T *>(src);
				return new (p)T(*srcDerived);
			}
			void destruct(void *p) const override {
				T *pDerived = static_cast<T *>(p);
				pDerived->~T();
			}
			static InstanceManager<TBase> *instance() {
				static InstanceManagerT<TBase, T> obj;
				return &obj;
			}
		};
	}


	template <class TBase, class... Ts>
	class StackBasedPolymophicValue {
	public:
		static_assert(sbpv_details::is_base_all_of<TBase, Ts...>::value, "Ts... must be inherits TBase");

		typedef StackBasedPolymophicValue<TBase, Ts...> Self;

		StackBasedPolymophicValue() {

		}
		~StackBasedPolymophicValue() {
			if (_manager) {
				_manager->destruct(p());
				_manager = nullptr;
			}
		}
		StackBasedPolymophicValue(const Self &rhs) {
			_manager = rhs._manager;
			if (_manager) {
				_manager->copyConstruct(p(), rhs.cp());
			}
		}
		StackBasedPolymophicValue(Self &&rhs) noexcept {
			_manager = rhs._manager;
			if (_manager) {
				_manager->copyConstruct(p(), rhs.cp());
			}
		}

		template <class T>
		StackBasedPolymophicValue(const T &value) {
			static_assert(sbpv_details::is_one_of<T, Ts...>::value, "T is must be any of Ts...");
			_manager = sbpv_details::InstanceManagerT<TBase, T>::instance();
			_manager->copyConstruct(p(), &value);
		}

		operator bool() const {
			return _manager != nullptr;
		}

		Self &operator=(std::nullptr_t) {
			if (_manager) {
				_manager->destruct(p());
				_manager = nullptr;
			}
			return *this;
		}

		template <class T>
		Self &operator=(const T &value) {
			static_assert(sbpv_details::is_one_of<T, Ts...>::value, "T is must be any of Ts...");

			if (_manager) {
				_manager->destruct(p());
			}
			_manager = sbpv_details::InstanceManagerT<TBase, T>::instance();
			_manager->copyConstruct(p(), &value);
			return *this;
		}
		Self &operator=(const Self &rhs) {
			if (this != &rhs) {
				if (_manager) {
					_manager->destruct(p());
				}
				_manager = rhs._manager;
				if (_manager) {
					_manager->copyConstruct(p(), rhs.cp());
				}
			}
			return *this;
		}
		TBase *operator->() {
			return _manager->ptr(p());
		}
		const TBase *operator->() const {
			return _manager->ptr(cp());
		}
	private:
		inline void *p() {
			return static_cast<void *>(&_storage);
		}
		inline const void *cp() const {
			return static_cast<const void *>(&_storage);
		}
	private:
		sbpv_details::InstanceManager<TBase> *_manager = nullptr;
		typename std::aligned_union<0, Ts...>::type _storage;
	};
}
