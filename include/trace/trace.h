//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <array>
#include <mutex>
#include <unordered_map>
#include <map>
#include <sstream>
#include <memory>
#include <vector>
#include <deque>
#include <algorithm>
#include <iomanip>
#include <thread>
#include <chrono>

#if defined(_WIN32)
#define NOMINMAX
#include <intrin.h>
#include <windows.h>
#endif

#if defined(__linux__)
#include <x86intrin.h>
#include <ctime>
#endif

#pragma once

namespace trace {

#if defined(__linux__)
    template< clockid_t id > struct clock_gettime_time_traits {
        using value_type = timespec;
        static value_type get() {
            value_type value;
            clock_gettime(id, &value);
            return value;
        }
        static double diff(const value_type& end, const value_type& begin) {
            uint64_t value = (end.tv_sec - begin.tv_sec) * 1e9 + (end.tv_nsec - begin.tv_nsec);
            return value;
        }
    };
#endif

#if defined(_WIN32)
    struct qpc_time_traits {
        using value_type = LARGE_INTEGER;
        static value_type get() {
            value_type value;
            QueryPerformanceCounter(&value);
            return value;
        }
        static double diff(const value_type& end, const value_type& begin) {
            value_type elapsed;
            elapsed.QuadPart = end.QuadPart - begin.QuadPart;
            elapsed.QuadPart *= 1e9;
            elapsed.QuadPart /= get_frequency().QuadPart;
            return elapsed.QuadPart;
        }

    private:
        static const value_type& get_frequency() {
            static value_type freq = []() {
                value_type value;
                QueryPerformanceFrequency(&value);
                return value;
            }();
            return freq;
        }
    };
#endif

    // This assumes constant inviarant TSC
    template< typename CTraits > struct rdtsc_time_traits {
        using value_type = uint64_t;

        static value_type get() {
            _mm_lfence();
            auto res = __rdtsc();
            _mm_lfence();
            return res;

            //unsigned dummy;
            //auto res = __rdtscp(&dummy);
            //_mm_lfence();
            //return res;
        }

        static double diff(const value_type& end, const value_type& begin) {
            return double(end - begin) / get_frequency();
        }

        static double get_frequency() {
            static double frequency = calculate_frequency();
            return frequency;
        }

        static uint64_t get_overhead() {
            static uint64_t overhead = calculate_overhead();
            return overhead;
        }
        
    private:
        static double calculate_frequency() {
            double freq = 0;
            const int N = 100;
            for (size_t i = 0; i < N; ++i) {
                auto cstart = CTraits::get();
                auto start = get();
	            for(size_t i = 0; i < 100000; ++i)
                    _mm_pause();
                auto stop = get();
                auto cstop = CTraits::get();
                freq += double(stop - start) / CTraits::diff(cstop, cstart);
            }

            return freq / N;
        }

        static uint64_t calculate_overhead() {
            const int N = 100;
            uint64_t x = 0, y = 0, min = -1;
      
            for (size_t i = 0; i < N; ++i) {
                x = get();
                y = get();
                min = std::min(min, y - x);
            }

            return min;
        }
    };

#if defined(_WIN32)
    using default_time_traits = rdtsc_time_traits< qpc_time_traits >;
#elif defined(__linux__)
    using default_time_traits = rdtsc_time_traits< clock_gettime_time_traits< CLOCK_MONOTONIC > >;
#else
    #error "don't know how to measure"
#endif

    struct frame {
        frame(const char* name) : _name(name) {}
        const char* name() const { return _name; }
    private:
        const char* _name;
    };

    struct frame_stack {
        struct hash {
            size_t operator()(const frame_stack& stack) const noexcept {
                size_t hval = 0;
                for (size_t i = 0; i < stack.addresses.size(); ++i)
                    hval ^= reinterpret_cast<uintptr_t>(stack.addresses[i]);
                return hval;
            }
        };

        bool operator == (const frame_stack& other) const { return addresses == other.addresses; }
        bool operator < (const frame_stack& other) const { return addresses < other.addresses; }

        std::deque< const frame* > addresses;
    };
    
    template< typename FrameData > struct frame_registry {
    private:
        using frame_records = std::unordered_map< frame_stack, FrameData, frame_stack::hash >;
        
        struct thread_data {
            thread_data(frame_records* records) : records(records) {}

            frame_stack stack;
            frame_records* records;
        };

        mutable std::mutex _mutex;
        std::vector< std::unique_ptr< frame_records > > _records;

        static frame_registry< FrameData > _instance;
        static thread_local thread_data _thread_data;
        static double _trace_overhead;

    public:
        static frame_registry< FrameData >& instance() { return _instance; }
            
        template< typename Fn > void for_each(Fn&& fn) const {
            auto records = merge();
            auto sorted_records = sort_tree(records);
            for_each(sorted_records, fn);
        }

        template< typename Fn > void for_each(const std::vector< FrameData* >& frames, Fn&& fn, size_t level = 0) const {
           for (auto& data: frames) {
                fn(*data, level);
                for_each(data->children(), fn, level + 1);
           }
        }

        FrameData* push_frame(const frame* frame) {
            auto& thread_data = _thread_data;
            thread_data.stack.addresses.push_back(frame);
            return get_frame_data(thread_data, thread_data.stack);
        }

        void pop_frame() {
            auto& thread_data = _thread_data;
            thread_data.stack.addresses.pop_back();
        }

        void clear() {
            for (auto& record : _records) {
                record->clear();
            }
        }

        double trace_overhead() { return _trace_overhead; }

    private:
        static std::string get_frame_name(const frame_stack& stack) {
            std::stringstream stream;
            stream << stack.addresses[0]->name();
            for(size_t i = 1; i < stack.addresses.size(); ++i) {
                stream << '/' << stack.addresses[i]->name();
            }
            return stream.str();
        }

        FrameData* get_frame_data(thread_data& data, const frame_stack& stack) const {
            FrameData* framedata = nullptr;
            auto& records = *data.records;
            auto it = records.find(stack);
            if (it != records.end()) {
                framedata = &it->second;
            } else {
                auto frame = stack.addresses.back();
                framedata = &(records.emplace(stack,
                    FrameData(get_parent_frame_data(data, stack), get_frame_name(stack), frame->name())).first->second);
            }
            return framedata;
        }

        FrameData* get_parent_frame_data(thread_data& data, const frame_stack& stack) const {
            FrameData* parent_data = nullptr;
            if (stack.addresses.size() > 1) {
                auto parent_stack(stack);
                parent_stack.addresses.pop_back();
                parent_data = get_frame_data(data, parent_stack);
            }
            return parent_data;
        }

        frame_records* create_frame_records() {
            std::lock_guard lock(_mutex);
            _records.emplace_back(new frame_records);
            return _records.back().get();
        }
        
        std::map< std::string, FrameData > merge() const {
            std::lock_guard lock(_mutex);
            std::map< std::string, FrameData > records;

            // Merge symbols together from different threads/callstacks
            for (auto& record : _records) {
                for (auto& [key, data] : *record) {
                    auto& mergeddata = records.emplace(data.long_name(), FrameData(data.parent(), data.long_name(), data.short_name())).first->second;
                    mergeddata.merge(data);
                }
            }

            // Create parent/child relationships according to merged symbols
            for (auto& [key, data] : records) {
                if (data.parent()) {
                    auto& parent = records[data.parent()->long_name()];
                    data.set_parent(&parent);
                    parent.add_child(&data);
                }
            }
            
            return records;
        }

        std::vector< FrameData* > sort_tree(std::map< std::string, FrameData >& records) const {
            auto compare = [&](const FrameData* lhs, const FrameData* rhs) { return lhs->time() > rhs->time(); };
            std::vector< FrameData* > result;
            for (auto& [key, data] : records) {
                data.sort_children(compare);
                if (!data.parent()) {
                    result.push_back(&data);
                }
            }

            std::sort(result.begin(), result.end(), compare);
            return result;
        }
    };

    template< typename FrameData > frame_registry<FrameData> frame_registry<FrameData>::_instance;
    template< typename FrameData > thread_local typename frame_registry<FrameData>::thread_data frame_registry<FrameData>::_thread_data(
        frame_registry<FrameData>::instance().create_frame_records());

    struct frame_data {
        frame_data() = default;

        frame_data(const frame_data* parent, std::string long_name, std::string short_name)
            : _long_name(long_name)
            , _short_name(short_name)
            , _parent(parent)
        {}

        void set_parent(const frame_data* parent) { _parent = parent; }
        const frame_data* parent() const { return _parent; }
        void add_child(frame_data* child) { _children.push_back(child); }
        const std::vector< frame_data* > children() const { return _children; }
        template< typename Fn > void sort_children(Fn&& fn) {
            std::sort(_children.begin(), _children.end(), fn);
            for (auto& child: _children) { child->sort_children(fn); }
        }

        const std::string& long_name() const { return _long_name; }
        const std::string& short_name() const { return _short_name; }

        void add_time(uint64_t value) { _time += value; _count += 1; }
        uint64_t time() const { return _time; }

        uint64_t exclusive_count() const { return _count; }

        uint64_t inclusive_count() const {
            uint64_t count = exclusive_count();
            for (auto& child : _children) { count += child->inclusive_count(); }
            return count;
        }

        void merge(const frame_data& other) {
            _time += other._time;
            _count += other._count;
        }
        
    private:
        // TODO: this does not belong here
        const frame_data* _parent {};
        std::vector< frame_data* > _children;
        std::string _long_name;
        std::string _short_name;

        uint64_t _time = 0;
        uint64_t _count = 0;
    };

    template< typename FrameData, typename TimeTraits = default_time_traits > struct frame_guard {
        frame_guard(const frame* frame): _frame_data(frame_registry< FrameData >::instance().push_frame(frame)) {
            _begin = TimeTraits::get();
        }

        ~frame_guard() {
            auto end = TimeTraits::get();
            _frame_data->add_time(TimeTraits::diff(end, _begin));
            frame_registry< FrameData >::instance().pop_frame();
        }

    private:
        FrameData* _frame_data;
        typename TimeTraits::value_type _begin;
    };

    template< typename FrameData > double frame_registry<FrameData>::_trace_overhead = []() {
        auto run = [](size_t n) {
            for (size_t i = 0; i < n; ++i) {
                static frame frame(__FUNCTION__);
                frame_guard< ::trace::frame_data > trace_guard(&frame);
            }
        };

        const size_t n = 100;
        run(n);

        auto begin = default_time_traits::get();
        run(n);
        auto end = default_time_traits::get();
        
        frame_registry<FrameData>::instance().clear();

        return default_time_traits::diff(end, begin) / n;
    }();
    
    struct stream_dumper {
        stream_dumper(std::ostream& stream, double scale = 1): _stream(stream), _scale(scale) {}

        void operator () (const frame_data& data, size_t level) {
            auto total_time = double(data.time()) * _scale / 1e6;
            auto avg_time = double(data.time()) / data.exclusive_count() / 1e6;
            auto overhead_count = std::max(uint64_t(0), data.inclusive_count() - data.exclusive_count());
            auto overhead_time = frame_registry< frame_data >::instance().trace_overhead() * overhead_count;
            auto cycles = uint64_t(data.time() * default_time_traits::get_frequency() / data.exclusive_count()) - default_time_traits::get_overhead();

            _stream << "CALLTRACE: " << indent(level) << std::fixed << std::setprecision(6)
                << data.short_name() << ": " << total_time
                << "ms (" << avg_time << "ms/call, count=" << data.exclusive_count();
            if (overhead_count) {
                _stream << ", err=" << std::setprecision(2) << (overhead_time * 100) / data.time() << "%";
            }
            if (cycles < 10000) {
                _stream << ", cycles=" << cycles;
            }
            _stream << ")" << std::endl;
        }

    private:
        std::string indent(size_t level) { return std::string(level * 4, ' '); }

        double _scale;
        std::ostream& _stream;
    };
}

#define _TRACE_CONCAT(a, b) _TRACE_CONCAT_(a, b)
#define _TRACE_CONCAT_(a, b) a ## b

#define _TRACE_NAME(name) _TRACE_CONCAT(name, __LINE__)

// TODO: 'static' here is important, as it allows us to use addresses of static variables
// instead of EIP to distinguish stacks. It also is very intrusive.
#define TRACE(...) \
  static ::trace::frame _TRACE_NAME(_entry_frame_)(__VA_ARGS__); \
  ::trace::frame_guard< ::trace::frame_data > _TRACE_NAME(_entry_)(&_TRACE_NAME(_entry_frame_))
