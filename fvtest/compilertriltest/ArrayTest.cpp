/*******************************************************************************
 * Copyright IBM Corp. and others 2023
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
 * [2] https://openjdk.org/legal/assembly-exception.html
 *
 * SPDX-License-Identifier: EPL-2.0 OR Apache-2.0 OR GPL-2.0 WITH Classpath-exception-2.0 OR LicenseRef-GPL-2.0 WITH Assembly-exception
 *******************************************************************************/

#include "OpCodeTest.hpp"
#include "default_compiler.hpp"
#include <vector>

static const int32_t returnValueForArraycmpGreaterThan = 2;
static const int32_t returnValueForArraycmpLessThan = 1;
static const int32_t returnValueForArraycmpEqual = 0;
/**
 * @brief TestFixture class for arraycmp test
 *
 * @details Used for arraycmp test with the arrays with same data.
 * The parameter is the length parameter for the arraycmp evaluator.
 */
class ArraycmpEqualTest : public TRTest::JitTest, public ::testing::WithParamInterface<int32_t> {};
/**
 * @brief TestFixture class for arraycmp test
 *
 * @details Used for arraycmp test which has mismatched element.
 * The first parameter is the length parameter for the arraycmp evaluator.
 * The second parameter is the offset of the mismatched element in the arrays.
 */
class ArraycmpNotEqualTest : public TRTest::JitTest, public ::testing::WithParamInterface<std::tuple<int32_t, int32_t>> {};

TEST_P(ArraycmpEqualTest, ArraycmpLenSameArray) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    /*
     * "address=0" parameter is needed for arraycmp opcode because "Call" property is set to the opcode.
     * We need "flags=15" parameter to set arrayCmpLen flag.
     * arrayCmpLen flag is defined as 0x8000, which is 1 << 15.
     */
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address] flags=[15]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(length, entry_point(&s1[0], &s1[0]));
}

TEST_P(ArraycmpEqualTest, ArraycmpLenEqualConstLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address] flags=[15]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(length, entry_point(&s1[0], &s2[0]));
}

TEST_P(ArraycmpEqualTest, ArraycmpLenEqualVariableLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address, Int32]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address] flags=[15]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iload parm=2)))))"
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *, int32_t)>();
    EXPECT_EQ(length, entry_point(&s1[0], &s2[0], length));
}

TEST_P(ArraycmpEqualTest, ArraycmpSameArray) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(returnValueForArraycmpEqual, entry_point(&s1[0], &s1[0]));
}

TEST_P(ArraycmpEqualTest, ArraycmpEqualConstLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(returnValueForArraycmpEqual, entry_point(&s1[0], &s2[0]));
}

TEST_P(ArraycmpEqualTest, ArraycmpEqualVariableLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = GetParam();
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address, Int32]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iload parm=2)))))"
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *, int32_t)>();
    EXPECT_EQ(returnValueForArraycmpEqual, entry_point(&s1[0], &s2[0], length));
}

INSTANTIATE_TEST_CASE_P(ArraycmpTest, ArraycmpEqualTest, ::testing::Range(1, 128));

TEST_P(ArraycmpNotEqualTest, ArraycmpLenNotEqualConstLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address] flags=[15]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x3f;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(offset, entry_point(&s1[0], &s2[0]));
}

TEST_P(ArraycmpNotEqualTest, ArraycmpLenNotEqualVariableLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address, Int32]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address] flags=[15]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iload parm=2)))))"
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x3f;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *, int32_t)>();
    EXPECT_EQ(offset, entry_point(&s1[0], &s2[0], length));
}

TEST_P(ArraycmpNotEqualTest, ArraycmpGreaterThanConstLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x81;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(returnValueForArraycmpGreaterThan, entry_point(&s1[0], &s2[0]));
}

TEST_P(ArraycmpNotEqualTest, ArraycmpGreaterThanVariableLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address, Int32]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iload parm=2)))))"
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x81;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *, int32_t)>();
    EXPECT_EQ(returnValueForArraycmpGreaterThan, entry_point(&s1[0], &s2[0], length));
}

TEST_P(ArraycmpNotEqualTest, ArraycmpLessThanConstLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iconst %d)))))",
      length
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x21;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *)>();
    EXPECT_EQ(returnValueForArraycmpLessThan, entry_point(&s1[0], &s2[0]));
}

TEST_P(ArraycmpNotEqualTest, ArraycmpLessThanVariableLen) {
    SKIP_ON_ARM(MissingImplementation);
    SKIP_ON_RISCV(MissingImplementation);

    auto length = std::get<0>(GetParam());
    auto offset = std::get<1>(GetParam());
    char inputTrees[1024] = {0};
    std::snprintf(inputTrees, sizeof(inputTrees),
      "(method return=Int32 args=[Address, Address, Int32]"
      "  (block"
      "    (ireturn"
      "      (arraycmp address=0 args=[Address, Address]"
      "        (aload parm=0)"
      "        (aload parm=1)"
      "        (iload parm=2)))))"
      );
    auto trees = parseString(inputTrees);

    ASSERT_NOTNULL(trees);

    Tril::DefaultCompiler compiler(trees);

    ASSERT_EQ(0, compiler.compile()) << "Compilation failed unexpectedly\n" << "Input trees: " << inputTrees;

    std::vector<unsigned char> s1(length, 0x5c);
    std::vector<unsigned char> s2(length, 0x5c);
    s1[offset] = 0x21;

    auto entry_point = compiler.getEntryPoint<int32_t (*)(unsigned char *, unsigned char *, int32_t)>();
    EXPECT_EQ(returnValueForArraycmpLessThan, entry_point(&s1[0], &s2[0], length));
}

static std::vector<std::tuple<int32_t, int32_t>> createArraycmpNotEqualParam() {
  std::vector<std::tuple<int32_t, int32_t>> v;
  /* Small arrays */
  for (int i = 1; i < 32; i++) {
    for (int j = 0; j < i; j++) {
      v.push_back(std::make_tuple(i, j));
    }
  }
  /* Variation of the offset of mismatched element in 128 bytes array */
  for (int i = 0; i < 128; i++) {
    v.push_back(std::make_tuple(128, i));
  }
  /* Medium size arrays with the mismatched element near the end of the arrays */
  for (int i = 120; i < 136; i++) {
    for (int j = 96; j < i; j++) {
      v.push_back(std::make_tuple(i, j));
    }
  }
  /* A large size array with the mismatched element near the end of the array */
  for (int i = 4000; i < 4096; i++) {
    v.push_back(std::make_tuple(4096, i));
  }
  return v;
}
INSTANTIATE_TEST_CASE_P(ArraycmpTest, ArraycmpNotEqualTest, ::testing::ValuesIn(createArraycmpNotEqualParam()));
