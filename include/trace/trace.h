//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

// TODO: 
// move sampling counter to the frame (not as easy as we will lose counts...)
// inject frame_registry
// get rid of time/count, use histogram (?)

#include <array>
#include <mutex>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <deque>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <thread>
#include <chrono>
#include <iostream>
#include <cmath>
#include <map>

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

#define TRACE_ENABLED

namespace trace {
    template< size_t N, size_t D = 1 > struct ratio {
        static_assert(D >= N);
        static constexpr double value() { return double(N) / D; }
    };

    using BeginRatio = ratio<0>;
    using EndRatio = ratio<3, 4>;

#if defined(__linux__)
    template< clockid_t id > struct clock_gettime_time_traits {
        using value_type = timespec;

        static value_type begin() { return get(); }
        static value_type end() { return get(); }

        static double diff(const value_type& end, const value_type& begin) {
            uint64_t value = (end.tv_sec - begin.tv_sec) * 1e9 + (end.tv_nsec - begin.tv_nsec);
            return value;
        }

    private:
        static value_type get() {
            value_type value;
            clock_gettime(id, &value);
            return value;
        }
    };
#endif

#if defined(_WIN32)
    struct qpc_time_traits {
        using value_type = LARGE_INTEGER;

        static value_type begin() { return get(); }
        static value_type end() { return get(); }
        
        static double diff(const value_type& end, const value_type& begin) {
            value_type elapsed;
            elapsed.QuadPart = end.QuadPart - begin.QuadPart;
            elapsed.QuadPart *= 1e9;
            elapsed.QuadPart /= get_frequency().QuadPart;
            return elapsed.QuadPart;
        }

    private:
        static value_type get() {
            value_type value;
            QueryPerformanceCounter(&value);
            return value;
        }

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
    // With mfence, the counter has fixed cost that can be substracted from TRACE().
    struct rdtscp_counter {
        using value_type = uint64_t;

        static value_type begin() { return get(); }
        static value_type end() { return get(); }

    private:
        static value_type get() {
            unsigned dummy;
            _mm_mfence();
            auto res = __rdtscp(&dummy);
            _mm_lfence();
            return res;
        }
    };

    struct rdtsc_counter {
        using value_type = uint64_t;

        static value_type begin() { return get(); }
        static value_type end() { return get(); }

    private:
        static value_type get() {
            _mm_mfence();
            _mm_lfence();
            auto res = __rdtsc();
            _mm_lfence();
            return res;
        }
    };

    template< typename CTraits > struct tsc_time_traits: rdtscp_counter {
        using value_type = uint64_t;
 
        static double diff(const value_type& end, const value_type& begin) {
            return double(end - begin) / get_frequency();
        }

        static double get_frequency() {
            static double frequency = calculate_frequency();
            return frequency;
        }

        static const std::tuple< uint64_t, double, uint64_t >& get_overhead() {
            static std::tuple< uint64_t, double, uint64_t > overhead = calculate_overhead();
            return overhead;
        }
        
    private:
        static double calculate_frequency() {
            const int N = 1024;
            std::array< double, N > data;
            for (size_t i = 0; i < N; ++i) {
                auto cstart = CTraits::begin();
                auto start = begin();
                for(size_t i = 0; i < 1024 * 32; ++i)
                    _mm_pause();
                auto stop = end();
                auto cstop = CTraits::end();
                data[i] = double(stop - start) / CTraits::diff(cstop, cstart);
            }

            std::sort(data.begin(), data.end(), [](auto lhs, auto rhs) { return lhs < rhs; });
            
            double total = 0;
            for (size_t i = BeginRatio::value() * N; i < EndRatio::value() * N; ++i) {
                total += data[i];
            }

            return total / ((EndRatio::value() - BeginRatio::value()) * N);
        }

        static std::tuple< uint64_t, double, uint64_t > calculate_overhead() {
            const int N = 1024*16;
            uint64_t x = 0, y = 0;

            uint64_t counter = 0;
            std::vector< uint64_t > data(N);
            for (size_t i = 0; i < N; ++i) {
                x = begin();
                y = end();
                data[i] = y - x;
            }

            std::sort(data.begin(), data.end(), [](uint64_t lhs, uint64_t rhs) { return lhs < rhs; });

            uint64_t total = 0;
            for (size_t i = BeginRatio::value() * N; i < EndRatio::value() * N; ++i) {
                total += data[i];
            }

            uint64_t min = data[BeginRatio::value() * N];
            uint64_t max = data[EndRatio::value() * N - 1];
            return { min, double(total) / ((EndRatio::value() - BeginRatio::value()) * N), max };
        }
    };

#if defined(_WIN32)
    using default_time_traits = tsc_time_traits< qpc_time_traits >;
#elif defined(__linux__)
    using default_time_traits = tsc_time_traits< clock_gettime_time_traits< CLOCK_MONOTONIC > >;
#else
    #error "don't know how to measure"
#endif
    
    // See: Accurate Measurement of Small Execution Times
    // https://uwaterloo.ca/embedded-software-group/sites/ca.embedded-software-group/files/uploads/files/ieee-esl-precise-measurements.pdf
    template < typename TimeTraits, typename Begin = BeginRatio, typename End = EndRatio, typename Fn > std::tuple< double, double > differential_measure(const Fn& fn) {
        // T1 = N1F + e
        // T2 = N2F + e
        // F = (T2-T1)/(N2-N1)
        // e = (N1T2-N2T1)/(N1-N2)

        const int N = 1024*64;
        std::vector< std::array< double, 2 > > data(N);

        int64_t begin, end, t1, t2;
        size_t n1 = 1, n2 = 10;
        
        for (size_t i = 0; i < N; ++i) {
        again:
            begin = TimeTraits::begin();
            fn();
            end = TimeTraits::end();
            t1 = end - begin;

            begin = TimeTraits::begin();
            for(size_t i = 0; i < n2; ++i)
                fn();
            end = TimeTraits::end();
            t2 = end - begin;
            if (t2 < t1)
                goto again;

            data[i][0] = double(t2 - t1) / (n2 - n1);
            data[i][1] = double(n1 * t2 - n2 * t1) / (n1 - n2);
        }

        std::sort(data.begin(), data.end(), 
            [](const std::array< double, 2 >& lhs, const std::array< double, 2 >& rhs) {
            return lhs[0] < rhs[0];
        });

        double e = 0, f = 0;
        for (size_t i = Begin::value() * N; i < End::value() * N; ++i) {
            f += data[i][0];
            e += data[i][1];
        }
        return { f / ((End::value() - Begin::value()) * N), e / ((End::value() - Begin::value()) * N) };
    }

    struct histogram {
        struct bin {
            uint64_t total;
            uint64_t count;

            void merge(const bin& other) {
                total += other.total;
                count += other.count;
            }
        };

        size_t index(uint64_t value) {
        #if defined(_WIN32)
            auto lz = __lzcnt64(value);
        #else
            auto lz =  __builtin_clzll(value);
        #endif
            return (sizeof(value) * 8 - 1) - lz;
        }

        void add(uint64_t value, size_t count){
            auto i = index(value);
            _data[i].total += value * count;
            _data[i].count += count;
        }

        void merge(const histogram& other) {
            for (size_t i = 0; i < _data.size(); ++i)
                _data[i].merge(other._data[i]);
        }

        std::vector< double > values() const {
            size_t n = 0;
            for (size_t i = 0; i < _data.size(); ++i) {
                n += _data[i].count;
            }
            std::vector< double > data(n);
            n = 0;
            for (size_t i = 0; i < _data.size(); ++i) {
                if (!_data[i].count)
                    continue;

                auto mean = double(_data[i].total) / _data[i].count;
                for (size_t j = 0; j < _data[i].count; ++j) {
                    data[n++] = mean;
                }
            }
            return data;
        }

        template< typename Begin = BeginRatio, typename End = EndRatio > double mean() const {
            auto data = values();
            double mean = 0;
            size_t count = 0;

            for (size_t i = std::floor(Begin::value() * data.size()); i < std::ceil(End::value() * data.size()); ++i) {
                mean += data[i];
                count++;
            }
            return count > 0 ? mean / count : 0;
        }

    private:
        std::array< bin, sizeof(uint64_t) * 8 > _data{};
    };

    inline double compensate(const std::tuple< uint64_t, double, uint64_t >& overhead, double value, size_t count) {
        // TODO: the compensation is crude
        const auto& [min, avg, max] = overhead;
        if (value > avg * count)
            return value - avg * count;
        if (value > min * count)
            return value - min * count;
        return 0;
    }

    struct frame {
        frame(const char* name) : _name(name) {}
        const char* name() const { return _name; }

    private:
        const char* _name;
    };

    template< typename FrameData > struct frame_tree;

    template< typename Value, size_t N > struct frame_tree_cache {
        Value* get(const frame* frame) {
            for (size_t i = 0; i < _frames.size(); ++i) {
                if (_frames[i] == frame) return _values[i];
            }
            return nullptr;
        }

        void put(const frame* frame, Value* value) {
            auto index = _index++ & (N - 1);
            _frames[index] = frame;
            _values[index] = value;
        }

        void clear() {
            for (size_t i = 0; i < _frames.size(); ++i) {
                _frames[i] = nullptr;
                _values[i] = nullptr;
            }
            _index = 0;
        }

        std::array< const frame*, N > _frames {};
        std::array< Value*, N > _values {};
        size_t _index {};
    };

    // This is 10% faster than std::deque<>
    template< typename T, size_t MaxSize > struct array_stack {
        void push_back(const T& value) { _values[_size++] = value; }
        void pop_back() { --_size; }
        const T& back() { return _values[_size - 1]; }
        size_t size() const { return _size; }
        const T& operator[](size_t n) const { return _values[n]; }

    private:
        std::array< T, MaxSize > _values;
        size_t _size{};
    };

    template< typename FrameData > struct frame_tree {
        frame_tree()
            : _frame()
            , _parent()
        {}

        frame_tree(const frame_tree< FrameData >* parent, const frame* frame, FrameData&& data)
            : _parent(parent)
            , _frame(frame)
            , _data(std::move(data))
        {}

        void clear() {
            _children.clear();
            _cache.clear();
        }

        const std::unordered_map< const frame*, std::unique_ptr< frame_tree< FrameData > > >& children() const { return _children; }
        FrameData* data() { return &_data; }
        const FrameData* data() const { return &_data; }

        template< size_t N > frame_tree< FrameData >* get(const frame* frame, const array_stack< frame_tree< FrameData >*, N >& stack) {
            frame_tree< FrameData >* frame_data = _cache.get(frame);
            if (frame_data)
                return frame_data;

            auto it = _children.find(frame);
            if (it == _children.end()) {
                it = _children.emplace(frame,
                    new frame_tree< FrameData >(this, frame,
                        FrameData(&_data, get_long_name(stack, frame->name()), frame->name()))).first;
            }

            frame_data = it->second.get();
            _cache.put(frame, frame_data);
            return frame_data;
        }

    private:
        template< size_t N > static std::string get_long_name(const array_stack< frame_tree< FrameData >*, N >& stack, const std::string& name) {
            std::stringstream stream;
            for (size_t i = 1; i < stack.size(); ++i) {
                stream << '/' << stack[i]->_frame->name();
            }
            if (stack.size() > 1)
                stream << '/';
            stream << name;
            return stream.str();
        }

        const frame* _frame;
        const frame_tree< FrameData >* _parent;
        frame_tree_cache< frame_tree< FrameData >, 8 > _cache;
        std::unordered_map< const frame*, std::unique_ptr< frame_tree< FrameData > > > _children;
        FrameData _data;
    };

    template< typename FrameData > struct frame_registry {
    private:
        struct thread_data {
            thread_data(frame_tree< FrameData >* root) {
                stack.push_back(root);
            }
            
            array_stack< frame_tree< FrameData >*, 1024 > stack;
        };

        mutable std::mutex _mutex;
        std::vector< std::unique_ptr< frame_tree< FrameData > > > _frames;

        static frame_registry< FrameData > _instance;
        static thread_local thread_data _thread_data;
        static double _trace_overhead;

    public:
        static frame_registry< FrameData >& instance() { return _instance; }
            
        template< typename Fn > void for_each(Fn&& fn) const {
            // TODO: this is really bad, sorted_records have pointers to records
            auto records = merge();
            auto sorted_records = sort(records);
            for_each(sorted_records[0]->children(), fn); // TODO
        }

        template< typename Fn > void for_each(const std::vector< FrameData* >& frames, Fn&& fn, size_t level = 0) const {
           for (auto& data: frames) {
                fn(*data, level);
                for_each(data->children(), fn, level + 1);
           }
        }

        FrameData* push_frame(const frame* frame) {
            auto tree = _thread_data.stack.back()->get(frame, _thread_data.stack);
            _thread_data.stack.push_back(tree);
            return tree->data();
        }

        void pop_frame() {
            _thread_data.stack.pop_back();
        }

        void clear() {
            for (auto& frame : _frames) {
                frame->clear();
            }
        }

        double trace_overhead() { return _trace_overhead; }

    private:
        frame_tree< FrameData >* create_frame() {
            std::lock_guard lock(_mutex);
            _frames.emplace_back(new frame_tree< FrameData >);
            return _frames.back().get();
        }

        void merge(std::map< std::string, FrameData >& records, const frame_tree< FrameData >& tree, size_t level) const {
            if (level > 0) {
                auto data = tree.data();
                auto& mergeddata = records.emplace(data->long_name(), FrameData(data->parent(), data->long_name(), data->short_name())).first->second;
                mergeddata.merge(*data);
            }

            for (auto& [_, tree] : tree.children()) {
                merge(records, *tree, level + 1);
            }
        }

        std::map< std::string, FrameData > merge() const {
            std::lock_guard lock(_mutex);
            std::map< std::string, FrameData > records;

            // Merge symbols together from different threads/callstacks
            for (auto& frame : _frames) {
                merge(records, *frame, 0);
            }

            // Create parent/child relationships according to merged symbols
            for (auto& [key, data] : records) {
                if (data.parent()) {
                    auto& parent = records[data.parent()->long_name()];
                    data.set_parent(&parent);
                    parent.add_child(&data);
                }
            }
            
            return records; // This expects the move to happen, is that guaranteed?
        }

        std::vector< FrameData* > sort(std::map< std::string, FrameData >& records) const {
            auto compare = [&](const FrameData* lhs, const FrameData* rhs) { return lhs->value() > rhs->value(); };
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
        frame_registry<FrameData>::instance().create_frame());

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

        void add_value(uint64_t value, uint64_t count) {
            _value += value * count;
            _count += count;
            _histogram.add(value, count);
        }
        uint64_t value() const { return _value; }

        uint64_t exclusive_count() const { return _count; }

        uint64_t inclusive_count() const {
            uint64_t count = exclusive_count();
            for (auto& child : _children) { count += child->inclusive_count(); }
            return count;
        }

        void merge(const frame_data& other) {
            _value += other._value;
            _count += other._count;
            _histogram.merge(other._histogram);
        }

        uint64_t increment_sampling_counter() {
            return _sampling_counter[0]++;
        }

        uint64_t reset_sampling_counter() {
            auto elapsed = _sampling_counter[0] - _sampling_counter[1];
            _sampling_counter[1] = _sampling_counter[0];
            return elapsed;
        }

    private:
        // TODO: this does not belong here
        const frame_data* _parent {};
        std::vector< frame_data* > _children;
        std::string _long_name;
        std::string _short_name;

        uint64_t _value = 0;
        uint64_t _count = 0;
        uint64_t _sampling_counter[2] = {0};
    public:
        histogram _histogram;
    };

    template< typename FrameData, typename TimeTraits = default_time_traits, size_t SamplingFreq = 32 > struct frame_guard {
        frame_guard(const frame* frame) {
            _frame_data = frame_registry< FrameData >::instance().push_frame(frame);
            _counter = _frame_data->increment_sampling_counter();
            if ((_counter & (SamplingFreq - 1)) == 0)
                _begin = TimeTraits::begin();
        }

        ~frame_guard() {
            if ((_counter & (SamplingFreq - 1)) == 0) {
                auto end = TimeTraits::end();
                if (end > _begin)
                    _frame_data->add_value(end - _begin, _frame_data->reset_sampling_counter());
            }
            frame_registry< FrameData >::instance().pop_frame();
        }

        static const std::tuple< uint64_t, double, uint64_t >& get_overhead() {
            static std::tuple< uint64_t, double, uint64_t > overhead = calculate_overhead();
            return overhead;
        }

    private:
        static std::tuple< uint64_t, double, uint64_t > calculate_overhead() {
            const int N = 1024*16;
            uint64_t x = 0, y = 0;

            uint64_t counter = 0;
            size_t n = 0;
            std::vector< uint64_t > data(N);
            for (size_t i = 0; i < N * SamplingFreq; ++i) {
                if ((counter & (SamplingFreq - 1)) == 0)
                    x = TimeTraits::begin();
                if ((counter & (SamplingFreq - 1)) == 0)
                    y = TimeTraits::end();
                if ((counter & (SamplingFreq - 1)) == 0) {
                    data[n++] = y - x;
                }
                counter++;
            }

            std::sort(data.begin(), data.end(), [](uint64_t lhs, uint64_t rhs) { return lhs < rhs; });

            uint64_t total = 0;
            for(size_t i = N / 4; i < 3 * N / 4; ++i) {
                total += data[i];
            }

            uint64_t min = data[N/4];
            uint64_t max = data[3*N/4];
            return { min, double(total) / (n / 2), max };
        }

        uint64_t _counter;
        FrameData* _frame_data;
        typename TimeTraits::value_type _begin {};
    };

    template< typename FrameData > double frame_registry<FrameData>::_trace_overhead = []() {
        auto run = [](){
            static frame frame(__FUNCTION__);
            frame_guard< ::trace::frame_data > trace_guard(&frame);
        };
        
        run();
        
        auto [f, e] = differential_measure< default_time_traits >(run);
        frame_registry<FrameData>::instance().clear();
        return f;
    }();
    
    struct stream_dumper {
        stream_dumper(std::ostream& stream, double scale = 1): _stream(stream), _scale(scale) {}

        void operator () (const frame_data& data, size_t level) {
            auto overhead = frame_guard< frame_data >::get_overhead();
            auto time_total = trace::compensate(overhead, data.value(), data.exclusive_count()) * _scale / default_time_traits::get_frequency() / 1e6;
            auto time_avg = time_total / data.exclusive_count();
            auto overhead_count = std::max(uint64_t(0), data.inclusive_count() - data.exclusive_count());
            auto overhead_cycles = frame_registry< frame_data >::instance().trace_overhead() * overhead_count;
            auto cycles_raw = data.value() / data.exclusive_count();
            auto cycles_mean = data._histogram.mean();
            auto cycles_compensated = trace::compensate(overhead, cycles_mean, 1);

            _stream << "CALLTRACE: " << indent(level) << std::fixed << std::setprecision(6)
                << data.short_name() << ": " << time_total
                << "ms (" << time_avg << "ms/call, count=" << data.exclusive_count();
            if (cycles_raw < 1000) {
                _stream << ", cycles(raw=" << std::setprecision(2) << cycles_raw << ", mean=" << cycles_mean << ", est=" << cycles_compensated << ")";
            }
            if (overhead_count) {
                _stream << ", err=" << std::setprecision(2) << (overhead_cycles * 100) / data.value() << "%";
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
#if defined(TRACE_ENABLED)
    #define TRACE(...) \
      static ::trace::frame _TRACE_NAME(_entry_frame_)(__VA_ARGS__); \
      ::trace::frame_guard< ::trace::frame_data > _TRACE_NAME(_entry_)(&_TRACE_NAME(_entry_frame_))
#else
    #define TRACE(...)
#endif
