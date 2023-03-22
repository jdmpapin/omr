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

   noRegMask         = 0x00000000,

   // GPR
   //
   eaxMask           = 0x00000001,
   ebxMask           = 0x00000002,
   ecxMask           = 0x00000004,
   edxMask           = 0x00000008,
   ediMask           = 0x00000010,
   esiMask           = 0x00000020,
   ebpMask           = 0x00000040,
   espMask           = 0x00000080,
   AvailableGPRMask  = 0x000000FF,

   // FPR
   //
   st0Mask           = 0x00000001,
   st1Mask           = 0x00000002,
   st2Mask           = 0x00000004,
   st3Mask           = 0x00000008,
   st4Mask           = 0x00000010,
   st5Mask           = 0x00000020,
   st6Mask           = 0x00000040,
   st7Mask           = 0x00000080,
   AvailableFPRMask  = 0x000000FF,

   // XMMR
   //
   xmm0Mask          = 0x00000001,
   xmm1Mask          = 0x00000002,
   xmm2Mask          = 0x00000004,
   xmm3Mask          = 0x00000008,
   xmm4Mask          = 0x00000010,
   xmm5Mask          = 0x00000020,
   xmm6Mask          = 0x00000040,
   xmm7Mask          = 0x00000080,
   AvailableXMMRMask = 0x000000FF,

   // Vector Mask Registers
   //
   k0Mask            = 0x00000001,
   k1Mask            = 0x00000002,
   k2Mask            = 0x00000004,
   k3Mask            = 0x00000008,
   k4Mask            = 0x00000010,
   k5Mask            = 0x00000020,
   k6Mask            = 0x00000040,
   k7Mask            = 0x00000080,
   AvailableKMask    = 0x000000FF,
