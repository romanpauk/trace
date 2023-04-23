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
#include <algorithm>
#include <iomanip>

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
    struct time_traits {
        using value_type = timespec;
        static void get(value_type& value) { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value); }
        static uint64_t diff(const value_type& end, const value_type& begin) {
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

    template< typename FrameData > struct frame_registry {
    private:
        static const int MaxDepth = 16;

        struct frame_stack {
            size_t size;
            std::array< const frame*, MaxDepth > addresses;
            
            bool operator == (const frame_stack& other) const {
                return addresses == other.addresses;
            }

            bool operator < (const frame_stack& other) const {
                return addresses < other.addresses;
            }
        };

        struct frame_stack_hash {
            size_t operator()(const frame_stack& stack) const noexcept {
                size_t hval = 0;
                for (size_t i = 0; i < MaxDepth; ++i)
                    hval ^= reinterpret_cast<uintptr_t>(stack.addresses[i]);
                return hval;
            }
        };

        using frame_records = std::unordered_map< frame_stack, FrameData, frame_stack_hash >;

        struct thread_data {
            thread_data(frame_records* records) : records(records) {}

            frame_stack stack{};
            frame_records* records;
        };

        mutable std::mutex _mutex;
        std::vector< std::unique_ptr< frame_records > > _records;

        static frame_registry< FrameData > _instance;
        static thread_local thread_data _thread_data;

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
            auto& data = _thread_data;
            auto level = data.stack.size++;
            data.stack.addresses[level] = frame;

            // Create frame data
            return get_frame_data(data, data.stack);
        }

        void pop_frame() {
            auto& data = _thread_data;
            data.stack.addresses[--data.stack.size] = 0;
        }

    private:
        static std::string get_frame_name(const std::array< const frame*, MaxDepth >& stack) {
            std::stringstream stream;
            stream << stack[0]->name();
            for(size_t i = 1; i < MaxDepth; ++i) {
                if (!stack[i])
                    break;
                stream << '/' << stack[i]->name();
            }
            return stream.str();
        }

        FrameData* get_frame_data(thread_data& data, const frame_stack& stack) const {
            frame_records& records = *data.records;
            auto level = stack.size - 1;
            FrameData* framedata = 0;
            if (framedata)
                return framedata;
            auto it = records.find(stack);
            if (it != records.end()) {
                framedata = &it->second;
            }
            else {
                auto frame = stack.addresses[stack.size - 1];
                framedata = &(records.emplace(stack,
                    FrameData(get_parent_frame_data(data, stack), get_frame_name(stack.addresses), frame->name())).first->second);
            }
            return framedata;
        }

        FrameData* get_parent_frame_data(thread_data& data, const frame_stack& stack) const {
            FrameData* parent_data = nullptr;
            if (stack.size > 1) {
                auto parent_stack(stack);
                parent_stack.addresses[--parent_stack.size] = 0;
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
                    auto& mergeddata = records.emplace(data.name(), FrameData(data.parent(), data.name(), data.short_name())).first->second;
                    mergeddata.merge(data);
                }
            }

            // Create parent/child relationships according to merged symbols
            for (auto& [key, data] : records) {
                if (data.parent()) {
                    auto& parent = records[data.parent()->name()];
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

        frame_data(const frame_data* parent, std::string name, std::string short_name)
            : _name(name)
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

        const std::string& name() const { return _name; }
        const std::string& short_name() const { return _short_name; }

        void add_time(uint64_t value) { _time += value; _count += 1; }
        uint64_t time() const { return _time; }
        uint64_t call_count() const { return _count; }

        void merge(const frame_data& other) {
            _time += other._time;
            _count += other._count;
        }
        
    private:
        // TODO: this does not belong here
        const frame_data* _parent {};
        std::vector< frame_data* > _children;
        std::string _name;
        std::string _short_name;

        uint64_t _time = 0;
        uint64_t _count = 0;
    };

    template< typename FrameData, typename TimeTraits = time_traits > struct frame_guard {
        frame_guard(const frame* frame): _frame_data(frame_registry< FrameData >::instance().push_frame(frame)) {
            TimeTraits::get(_begin);
        }

        ~frame_guard() {
            typename TimeTraits::value_type end;
            TimeTraits::get(end);
            _frame_data->add_time(TimeTraits::diff(end, _begin));
            frame_registry< FrameData >::instance().pop_frame();
        }

    private:
        FrameData* _frame_data;
        typename TimeTraits::value_type _begin;
    };

    struct stream_dumper {
        stream_dumper(std::ostream& stream, size_t div = 1): _stream(stream), _div(div) {}

        void operator () (const frame_data& data, size_t level) {
            _stream << "CALLTRACE: " << indent(level) << data.short_name() << ": " << std::fixed << std::setprecision(6) << double(data.time()) / 1e6 << "ms (" << double(data.time())/1.e6/data.call_count() << "ms/call), count " << data.call_count() << std::endl;
        }

    private:
        std::string indent(size_t level) { return std::string(level * 4, ' '); }

        size_t _div;
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
  
