/* jbwlib/perf.h
 *   this file is part of jbwlib, a personal C routine library
 *
 * Copyright (C) 2014 Justin Willmert <justin@jdjlab.com>
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Configuration options:
 *
 * DEBUG_NOFUNC
 *     If defined, then log messages will not include the function name in
 *     the message output.
 *
 * DEBUG_NOCOLOR
 *     If defined, then terminal color control characters will not be emitted.
 */

#ifndef LIBJMERT_DEBUG_H
#define LIBJMERT_DEBUG_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
    #include <stdio.h>

    #define GCC_VERSION (__GNUC__ * 10000 \
                         + __GNUC_MINOR__* 100 \
                         + __GNUC_PATCHLEVEL__)

    // Get the correct function name to be used
    #if GCC_VERSION > 30200
        #define DBGLOG_FUNC __PRETTY_FUNCTION__
    #else
        #define DBGLOG_FUNC __func__
    #endif

    #ifdef DEBUG_NOCOLOR
        #define COLOR_BEGIN ""
        #define COLOR_END   ""
    #else /* DEBUG_NOCOLOR */
        #define COLOR_BEGIN "\x1b[31m"
        #define COLOR_END   "\x1b[0m"
    #endif

    // If function annotations are disabled,
    #ifdef DEBUG_NOFUNC

        #define xdbglog(cb, ce, msg, file, line, ...) \
            fprintf(stderr, cb "(%s:%d): " msg ce, \
                    file, line, ##__VA_ARGS__)
        #define dbglog(msg, ...) \
            xdbglog(COLOR_BEGIN, COLOR_END, msg, __FILE__, \
                    __LINE__, ##__VA_ARGS__)

    #else /* DEBUG_NOFUNC */

        #define xdbglog(cb, ce, msg, func, file, line, ...) \
            fprintf(stderr, cb "%s (%s:%d): " msg ce, \
                    func, file, line, ##__VA_ARGS__)
        #define dbglog(msg, ...) \
            xdbglog(COLOR_BEGIN, COLOR_END, msg, DBGLOG_FUNC, __FILE__, \
                    __LINE__, ##__VA_ARGS__)

    #endif

#else /* NDEBUG */

    /* Give an empty definition of all macros */
    #define dbglog(msg, ...)

#endif /* NDEBUG */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBJMERT_DEBUG_H */
