//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <gtest/gtest.h>

#include <trace/trace.h>

#include <iostream>

TEST(trace_test, test) {
    {
        {
            for(int i = 0; i < 10; ++i) {
                TRACE("a");
            }
        }
        TRACE("b");
        TRACE("c");
    }

    trace::registry< trace::callstack_storage >::instance().for_each([](const trace::callstack_storage& data) {
        std::cout << data.name() << ", " << data.total() / 1e6 / data.count() << "ms/call" << std::endl;
    });
}
