/*******************************************************************************
 * Copyright IBM Corp. and others 1991
 *
 * This program and the accompanying materials are made available under
 * the terms of the Eclipse Public License 2.0 which accompanies this
 * distribution and is available at http://eclipse.org/legal/epl-2.0
 * or the Apache License, Version 2.0 which accompanies this distribution
 * and is available at https://www.apache.org/licenses/LICENSE-2.0.
 *
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the
 * Eclipse Public License, v. 2.0 are satisfied: GNU General Public License,
 * version 2 with the GNU Classpath Exception [1] and GNU General Public
 * License, version 2 with the OpenJDK Assembly Exception [2].
 *
 * [1] https://www.gnu.org/software/classpath/license.html
 * [2] http://openjdk.java.net/legal/assembly-exception.html
 *
 * SPDX-License-Identifier: EPL-2.0 OR Apache-2.0 OR GPL-2.0 WITH Classpath-exception-2.0 OR LicenseRef-GPL-2.0 WITH Assembly-exception
 *******************************************************************************/

/**
 * @file
 * @ingroup PortTest
 * @brief Verify TTY operations
 *
 * Exercise the omrtty_daemonize API for port library functions in @ref omrtty.c
 *
 *
 */
#include <string.h>
#include <stdio.h>
#if defined(OMR_OS_WINDOWS)
#include <windows.h>
#endif /* defined(OMR_OS_WINDOWS) */

#include "testHelpers.hpp"
#include "omrport.h"

/**
 * Verify that a call to omrtty_printf/omrtty_vprintf after a call to omrtty_daemonize does not crash.
 *
 * Call omrtty_daemonize, then call omrtty_printf/omrtty_vprintf. If we crash or print anything, then the test has failed.
 */
#if defined(OMR_OS_WINDOWS)
TEST(PortTtyTest, tty_daemonize_test)
{
	OMRPORT_ACCESS_FROM_OMRPORT(portTestEnv->getPortLibrary());
	const char *format = "Test string %d %d %d\n";


	/* Get current output handles. */
	HANDLE stdoutOriginal = GetStdHandle(STD_OUTPUT_HANDLE);
	ASSERT_NE(stdoutOriginal, INVALID_HANDLE_VALUE);
	HANDLE stderrOriginal = GetStdHandle(STD_ERROR_HANDLE);
	ASSERT_NE(stderrOriginal, INVALID_HANDLE_VALUE);


	/* Redirect output to a test handle and restart TTY to get the port functions to use it. */
	HANDLE stdoutWrite;
	HANDLE stdoutRead;
	ASSERT_TRUE(CreatePipe(&stdoutRead, &stdoutWrite, NULL, 0));
	SetStdHandle(STD_OUTPUT_HANDLE, stdoutWrite);
	SetStdHandle(STD_ERROR_HANDLE, stdoutWrite);
	omrtty_shutdown();
	omrtty_startup();
	/* restore the original handlers as soon as possible to prevent gtest framework from attempting to use them */
	SetStdHandle(STD_OUTPUT_HANDLE, stdoutOriginal);
	SetStdHandle(STD_ERROR_HANDLE, stderrOriginal);

	omrtty_printf("pre-daemonize\n");
	omrtty_daemonize();
	/* outputComment uses omrtty_vprintf */
	outputComment(OMRPORTLIB, format, 1, 2, 3);
	omrtty_printf("post-daemonize\n");

	char readBuffer[EsMaxPath] = {0};
	DWORD bytesRead;
	EXPECT_TRUE(ReadFile(stdoutRead, readBuffer, EsMaxPath - 1, &bytesRead, NULL));

	/* Restart TTY to resume output for remaining tests. */
	omrtty_shutdown();
	omrtty_startup();

	/* Check the output. */
	outputComment(OMRPORTLIB, "Captured test output: %s\n", readBuffer);


	EXPECT_NE(strstr(readBuffer, "pre-daemonize"), (char*)NULL) << "'pre-daemonize' not found in output";
	EXPECT_EQ(strstr(readBuffer, "post-daemonize"), (char*)NULL) << "Unexpected 'post-daemonize' found in output";
	EXPECT_EQ(strstr(readBuffer, "Test string"), (char*)NULL) << "Unexpected 'Test string' found in output";

}
#endif /* defined(OMR_OS_WINDOWS) */
