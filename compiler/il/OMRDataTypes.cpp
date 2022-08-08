/*******************************************************************************
 * Copyright (c) 2000, 2022 IBM Corp. and others
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

#include "il/DataTypes.hpp"

#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "env/jittypes.h"
#include "env/TRMemory.hpp"
#include "il/DataTypes.hpp"
#include "il/ILOps.hpp"
#include "infra/Assert.hpp"
#include "infra/String.hpp"

// When adding new types also update pDataTypeNames[] in ras/Tree.cpp

bool
OMR::DataType::canGetMaxPrecisionFromType()
   {
   switch (self()->getDataType())
      {
      case TR::Int8:
      case TR::Int16:
      case TR::Int32:
      case TR::Int64:
         return true;
      default:
         return false;
      }
   }

int32_t
OMR::DataType::getMaxPrecisionFromType()
   {
   switch (self()->getDataType())
      {
      case TR::Int8: return TR::getMaxSignedPrecision<TR::Int8>();
      case TR::Int16: return TR::getMaxSignedPrecision<TR::Int16>();
      case TR::Int32: return TR::getMaxSignedPrecision<TR::Int32>();
      case TR::Int64: return TR::getMaxSignedPrecision<TR::Int64>();
      default:
         TR_ASSERT(false, "Unsupported data type in getMaxPrecisionFromType\n");
         return 0;
      }
   }

TR::DataType
OMR::DataType::getVectorIntegralType()
   {
   if (!self()->isVector()) return TR::NoType;

   switch(self()->getVectorElementType())
      {
      case TR::Int8:
      case TR::Int16:
      case TR::Int32:
      case TR::Int64: return self()->getDataType();
      case TR::Float: return TR::DataType::createVectorType(TR::Int32, self()->getVectorLength());
      case TR::Double: return TR::DataType::createVectorType(TR::Int64, self()->getVectorLength());
      default:
         return TR::NoType;
         break;
      }
   }

TR::DataType
OMR::DataType::vectorToScalar()
   {
   return self()->getVectorElementType();
   }

TR::DataType
OMR::DataType::scalarToVector(TR::VectorLength length)
   {
   TR::DataTypes type = self()->getDataType();

   return createVectorType(type, length);
   }

TR::DataType
OMR::DataType::getIntegralTypeFromPrecision(int32_t precision)
   {
   if (precision < 1 || precision >= TR::getMaxSignedPrecision<TR::Int64>())
      return TR::NoType;
   else if (precision < TR::getMaxSignedPrecision<TR::Int8>())
      return TR::Int8;
   else if (precision < TR::getMaxSignedPrecision<TR::Int16>())
      return TR::Int16;
   else if (precision < TR::getMaxSignedPrecision<TR::Int32>())
      return TR::Int32;
   else
      return TR::Int64;
   }

TR::DataType
OMR::DataType::getFloatTypeFromSize(int32_t size)
   {
   TR::DataType type = TR::NoType;
   switch (size)
      {
      case 4:  type = TR::Float; break;
      case 8:  type = TR::Double; break;
      default: TR_ASSERT(false,"unexpected size %d for float type\n",size);
      }
   return type;
   }


static int32_t OMRDataTypeSizes[] =
   {
   0,                 // TR::NoType
   1,                 // TR::Int8
   2,                 // TR::Int16
   4,                 // TR::Int32
   8,                 // TR::Int64
   4,                 // TR::Float
   8,                 // TR::Double
   sizeof(intptr_t), // TR::Address
   0,                 // TR::Aggregate
   };

static_assert(TR::NumOMRTypes == (sizeof(OMRDataTypeSizes) / sizeof(OMRDataTypeSizes[0])), "OMRDataTypeSizes is not the correct size");

int32_t
OMR::DataType::getSize(TR::DataType dt)
   {
   if (dt.isVector())
      {
      switch (dt.getVectorLength())
         {
         case TR::VectorLength64:  return 8;
         case TR::VectorLength128: return 16;
         case TR::VectorLength256: return 32;
         case TR::VectorLength512: return 64;
         default:
            TR_ASSERT_FATAL(false, "Incorrect Vector Length\n");
         }
      }

   TR_ASSERT(dt < TR::NumOMRTypes, "dataTypeSizeMap called on unrecognized data type");
   return OMRDataTypeSizes[dt];
   }

void
OMR::DataType::setSize(TR::DataType dt, int32_t newSize)
   {
   if (dt.isVector()) return;

   TR_ASSERT(dt < TR::NumOMRTypes, "setDataTypeSizeInMap called on unrecognized data type");
   OMRDataTypeSizes[dt] = newSize;
   }


static const char * OMRDataTypeNames[TR::NumAllTypes] =
   {
   "NoType",
   "Int8",
   "Int16",
   "Int32",
   "Int64",
   "Float",
   "Double",
   "Address",
   "Aggregate",
   // vector names will be created at runtime
   };

#define MAX_TYPE_NAME_LENGTH 20

bool
OMR::DataType::initVectorNames()
   {
   int32_t numVectorTypes = TR::NumAllTypes - TR::NumScalarTypes;
   char *names = (char*)TR_Memory::jitPersistentAlloc(MAX_TYPE_NAME_LENGTH*numVectorTypes*sizeof(char));
   char *name = names;

   for (int32_t i = TR::NumScalarTypes; i < TR::NumAllTypes; i++)
      {
      TR::DataType dt((TR::DataTypes)i);
      TR_ASSERT_FATAL(dt.isVector(), "Should be a vector type");
      TR::snprintfNoTrunc(name, MAX_TYPE_NAME_LENGTH, "Vector%d%s", getSize(dt)*8, getName(dt.getVectorElementType()));
      OMRDataTypeNames[dt] = name;
      name += MAX_TYPE_NAME_LENGTH;
      }
   return true;
   }

const char *
OMR::DataType::getName(TR::DataType dt)
   {
   if (dt.isVector())
      {
      // to avoid any race conditions, initialize all vector names once,
      // as soon as first one is requested
      static bool staticallyInitialized = initVectorNames();
      TR_ASSERT_FATAL(staticallyInitialized && (OMRDataTypeNames[dt] != NULL), "Vector names should've been initialized");
      }

   TR_ASSERT(dt.isOMRDataType(), "Name requested for a non-OMR datatype\n");
   return OMRDataTypeNames[dt];
   }

TR::DataType
OMR::DataType::getTypeFromName(const char *name)
   {
   for (int32_t i = 1 ; i < TR::NumAllTypes; i++)
      {
      TR::DataType dt = (TR::DataTypes)i;

      if (!dt.isOMRDataType()) continue;

      if (strcmp(name, getName(dt)) == 0)
         return  dt;
      }

   return TR::NoType;
   }

const char *
OMR::DataType::toString() const
   {
   return TR::DataType::getName(self()->getDataType());
   }

void FloatingPointLimits::setMaxFloat()
   {
   int32_t f = FLOAT_POS_INFINITY;
   _maxFloat = 0.0f;
   memcpy(&_maxFloat, &f, sizeof(f));
   }

void  FloatingPointLimits::setMinFloat()
   {
   int32_t f = FLOAT_NEG_INFINITY;
   _minFloat = 0.0f;
   memcpy(&_minFloat, &f, sizeof(f));
   }
void  FloatingPointLimits::setMaxDouble()
   {
   uint64_t d = DOUBLE_POS_INFINITY;
   _maxDouble = 0.0;
   memcpy(&_maxDouble, &d, sizeof(d));
   }

void  FloatingPointLimits::setMinDouble()
   {
   int64_t d = DOUBLE_NEG_INFINITY;
   _minDouble = 0.0;
   memcpy(&_minDouble, &d, sizeof(d));
   }

namespace TR{
   static FloatingPointLimits fpLimits;

   float getMaxFloat() { return fpLimits.getMaxFloat(); }
   float getMinFloat() { return fpLimits.getMinFloat(); }
   double getMaxDouble() { return fpLimits.getMaxDouble(); }
   double getMinDouble() { return fpLimits.getMinDouble(); }
}
