/*******************************************************************************
 * Copyright IBM Corp. and others 2019
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

#include "env/CPU.hpp"
#include "env/CompilerEnv.hpp"
#include "env/jittypes.h"
#include "omrport.h"

bool
OMR::ARM64::CPU::isTargetWithinUnconditionalBranchImmediateRange(intptr_t targetAddress, intptr_t sourceAddress)
   {
   intptr_t range = targetAddress - sourceAddress;
   return range <= self()->maxUnconditionalBranchImmediateForwardOffset() &&
          range >= self()->maxUnconditionalBranchImmediateBackwardOffset();
   }

bool
OMR::ARM64::CPU::supportsFeature(uint32_t feature)
   {
   if (TR::Compiler->omrPortLib == NULL)
      {
      return false;
      }

   OMRPORT_ACCESS_FROM_OMRPORT(TR::Compiler->omrPortLib);
   return (TRUE == omrsysinfo_processor_has_feature(&_processorDescription, feature));
   }

const char*
OMR::ARM64::CPU::getProcessorName()
   {
   const char* returnString = "";
   switch(_processorDescription.processor)
      {
      case OMR_PROCESSOR_ARM64_UNKNOWN:
         returnString = "Unknown ARM64 processor";
         break;
      case OMR_PROCESSOR_ARM64_V8_A:
         returnString = "ARMv8-A processor";
         break;
      default:
         returnString = "Unknown ARM64 processor";
         break;
      }
   return returnString;
   }
