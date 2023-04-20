//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <atomic>
#include <mutex>
#include <deque>

#if defined(_WIN32)
#include <windows.h>
#endif

#if defined(__linux__)
#include <ctime>
#endif

#pragma once

namespace trace {
#if defined(_WIN32)
    struct TimeTraits {
        using value_type = LARGE_INTEGER;

        static void get(value_type& value) { QueryPerformanceCounter(&value); }
        
        static uint64_t diff(const value_type& end, const value_type& begin) {
            value_type elapsed;
            elapsed.QuadPart = end.QuadPart - begin.QuadPart;
            elapsed.QuadPart *= 1e9;
            elapsed.QuadPart /= frequency().QuadPart;
            return elapsed.QuadPart;
        }

    private:
        static const value_type& frequency() {
            static value_type freq = []()
            {
                value_type value;
                QueryPerformanceFrequency(&value);
                return value;
            }();
            return freq;
        }
    };
#else
    struct TimeTraits {
        using value_type = timespec;

        static void get(value_type& value) { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value); }

        static uint64_t diff(const value_type& end, const value_type& begin) {
            uint64_t value = (end.tv_sec - begin.tv_sec) * 1e9 + (end.tv_nsec - begin.tv_nsec);
            return value;
        }
    };
#endif

    template< typename Storage > struct registry {
        static registry< Storage >& instance() {
            static registry< Storage > storage;
            return storage;
        }

        void add(Storage* storage){
            std::lock_guard lock(_mutex);
            _storages.push_back(storage);
        }

        template< typename Fn > void for_each(Fn&& fn) {
            for (auto& storage : _storages) {
                fn(*storage);
            }
        }

    private:
        std::mutex _mutex;
        std::deque< Storage* > _storages;
    };

    struct callstack_storage {
        callstack_storage(const char* name): _name(name) {
            registry< callstack_storage >::instance().add(this);
        }

        const char* name() const { return _name; }

        void add(uint64_t value) { _total += value; _count += 1; }
        uint64_t total() const { return _total; }
        uint64_t count() const { return _count; }

    private:
        const char* _name;
        std::atomic< uint64_t > _total = 0;
        std::atomic< uint64_t > _count = 0;
    };

    template< typename TimeTraits = TimeTraits > struct callstack_entry {
        callstack_entry(callstack_storage* storage)
            : _storage(storage) {
            TimeTraits::get(_begin);
        }

        ~callstack_entry() {
            typename TimeTraits::value_type end;
            TimeTraits::get(end);
            
            _storage->add(TimeTraits::diff(end, _begin));
        }

    private:
        typename TimeTraits::value_type _begin;
        callstack_storage* _storage;
    };
}

#define _TRACE_CONCAT(a, b) _TRACE_CONCAT_(a, b)
#define _TRACE_CONCAT_(a, b) a ## b

#define _TRACE_NAME(name) _TRACE_CONCAT(name, __LINE__)

#define TRACE(...) \
  static ::trace::callstack_storage _TRACE_NAME(_trace_storage_)(__VA_ARGS__); \
  ::trace::callstack_entry<> _TRACE_NAME(_trace_entry_)(&_TRACE_NAME(_trace_storage_))
  
