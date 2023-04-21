//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <array>
#include <mutex>
#include <unordered_map>
#include <sstream>
#include <memory>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

#if defined(__linux__)
#include <ctime>
#endif

#pragma once

namespace trace {

#if defined(__linux__)
    struct time_traits
    {
        using value_type = timespec;

        static void get(value_type& value) { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value); }

        static uint64_t diff(const value_type& end, const value_type& begin)
        {
            uint64_t value = (end.tv_sec - begin.tv_sec) * 1e9 + (end.tv_nsec - begin.tv_nsec);
            return value;
        }
    };
#endif

#if defined(_WIN32)
    struct time_traits {
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
#endif

    struct frame {
        frame(const char* name) : _name(name) {}
        const char* name() const { return _name; }
    private:
        const char* _name;
    };

    template< typename FrameData, size_t MaxDepth, size_t N > struct frame_data_cache {
        struct FrameCache {
            std::array< const frame*, N > frames;
            std::array< FrameData*, N > data;
            int index;
        };

        FrameData* get(size_t level, const frame* frame) const {
            auto& cache = _cache[level];
            for (size_t i = 0; i < N; ++i) {
                if (cache.frames[i] == frame)
                    return cache.data[i];
            }
            return nullptr;
        }

        void put(size_t level, const frame* frame, FrameData* data) {
            auto& cache = _cache[level];
            cache.frames[cache.index] = frame;
            cache.data[cache.index] = data;
            cache.index = (cache.index + 1) & (N - 1);
        }

    private:
        std::array< FrameCache, MaxDepth > _cache {};
    };

    template< typename FrameData > struct frame_registry {
    private:
        static const int MaxDepth = 8;

        struct frame_stack_hash {
            size_t operator()(const std::array< const frame*, MaxDepth >& value) const noexcept {
                size_t hval = 0;
                for (size_t i = 0; i < MaxDepth; ++i)
                    hval ^= reinterpret_cast<uintptr_t>(value[i]);
                return hval;
            }
        };

        using frame_records = std::unordered_map< std::array< const frame*, MaxDepth >, FrameData, frame_stack_hash >;

    public:
        static frame_registry< FrameData >& instance() { return _instance; }
            
        template< typename Fn > void for_each(Fn&& fn) const {
            auto records = merge();
            for (auto& [key, data]: records) {
                fn(data);
            }
        }

        void push_frame(const frame* frame) {
            auto& data = _thread_data;
            data.stack.addresses[data.stack.size++] = frame;
        }

        FrameData& pop_frame() {
            auto& data = _thread_data;
            auto frame = data.stack.addresses[--data.stack.size];
            FrameData* framedata = data.cache.get(data.stack.size, frame);
            if (!framedata) {
                framedata = get_frame_data(*data.records, data.stack.addresses);
                _thread_data.cache.put(data.stack.size, frame, framedata);
            }
            data.stack.addresses[data.stack.size] = 0;
            return *framedata;
        }

    private:
        static std::string get_frame_name(const std::array< const frame*, MaxDepth >& stack) {
            std::stringstream stream;
            stream << stack[0]->name();
            for(size_t i = 1; i < MaxDepth && stack[i]; ++i)
                stream << '/' << stack[i]->name();
            return stream.str();
        }

        FrameData* get_frame_data(frame_records& records, const std::array< const frame*, MaxDepth >& stack) {
            auto it = records.find(stack);
            if (it != records.end()) {
                return &it->second;
            }
            else {
                return &(records.emplace(stack, FrameData(get_frame_name(stack))).first->second);
            }
        }

        frame_records* create_frame_records() {
            std::lock_guard lock(_mutex);
            _records.emplace_back(new frame_records);
            return _records.back().get();
        }
        
        frame_records merge() const {
            std::lock_guard lock(_mutex);
            frame_records records;
            for (auto& record : _records) {
                for (auto& [key, data] : *record) {
                    auto& mergeddata = records.emplace(key, data.name()).first->second;
                    mergeddata.merge(data);
                }
            }
            return records;
        }

        struct frame_stack {
            size_t size;
            std::array< const frame*, MaxDepth > addresses;
        };

        struct thread_data {
            thread_data(frame_records* records): records(records) {}

            frame_stack stack {};
            frame_data_cache< FrameData, MaxDepth, 4 > cache;
            frame_records* records;
        };

        mutable std::mutex _mutex;
        std::vector< std::unique_ptr< frame_records > > _records;

        static frame_registry< FrameData > _instance;
        static thread_local thread_data _thread_data;
    };

    template< typename FrameData > frame_registry<FrameData> frame_registry<FrameData>::_instance;
    template< typename FrameData > thread_local typename frame_registry<FrameData>::thread_data frame_registry<FrameData>::_thread_data(
        frame_registry<FrameData>::instance().create_frame_records());

    struct frame_data {
        frame_data(std::string name): _name(name) {}

        const std::string& name() const { return _name; }

        void add_time(uint64_t value) { _time += value; _count += 1; }
        uint64_t time() const { return _time; }
        uint64_t call_count() const { return _count; }

        void merge(const frame_data& other) {
            _time += other._time;
            _count += other._count;
        }
        
    private:
        std::string _name;
        uint64_t _time = 0;
        uint64_t _count = 0;
    };

    template< typename FrameData, typename TimeTraits = time_traits > struct frame_guard {
        frame_guard(const frame* frame) {
            frame_registry< FrameData >::instance().push_frame(frame);
            TimeTraits::get(_begin);
        }

        ~frame_guard() {
            typename TimeTraits::value_type end;
            TimeTraits::get(end);
            auto& frame_data = frame_registry< FrameData >::instance().pop_frame();
            frame_data.add_time(TimeTraits::diff(end, _begin));
        }

    private:
        typename TimeTraits::value_type _begin;
    };

    struct stream_dumper {
        stream_dumper(std::ostream& stream): _stream(stream) {}

        void operator () (const frame_data& data) {
            _stream << data.name() << ", " << data.time() / 1e6 / data.call_count() << "ms/call" << std::endl;
        }

    private:
        std::ostream& _stream;
    };
}

#define _TRACE_CONCAT(a, b) _TRACE_CONCAT_(a, b)
#define _TRACE_CONCAT_(a, b) a ## b

#define _TRACE_NAME(name) _TRACE_CONCAT(name, __LINE__)

#define TRACE(...) \
  ::trace::frame _TRACE_NAME(_entry_frame_)(__VA_ARGS__); \
  ::trace::frame_guard< ::trace::frame_data > _TRACE_NAME(_entry_)(&_TRACE_NAME(_entry_frame_))
  
