/*******************************************************************************
 * Copyright (c) 2021, 2021 IBM Corp. and others
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

#ifndef IDT_INCL
#define IDT_INCL

#include "compile/Compilation.hpp"
#include "optimizer/CallInfo.hpp"
#include "optimizer/abstractinterpreter/IDTNode.hpp"
#include "env/Region.hpp"
#include "env/VerboseLog.hpp"
#include "il/ResolvedMethodSymbol.hpp"
#include <queue>

namespace TR {

/**
 * IDT stands for Inlining Dependency Tree
 * It is a structure that holds all candidate methods to be inlined.
 *
 * The parent-child relationship in the IDT corresponds to the caller-callee relationship.
 */
class IDT
   {
   public:
   IDT(TR::Region& region, TR_CallTarget*, TR::ResolvedMethodSymbol* symbol, uint32_t budget, TR::Compilation* comp);

   TR::IDTNode* getRoot() { return _root; }

   TR::Region& getRegion() { return _region; }

   void addCost(uint32_t cost) { _totalCost += cost; }
   uint32_t getTotalCost() { return _totalCost; }

   /**
    * @brief Get the total number of nodes in this IDT.
    *
    * @return the total number of node
    */
   uint32_t getNumNodes() { return _nextIdx + 1; }

   /**
    * @brief Get the next avaible IDTNode index.
    *
    * @return the next index
    */
   int32_t getNextGlobalIDTNodeIndex() { return _nextIdx; }

   /**
    * @brief Increase the next available IDTNode index by 1.
    * This should only be called when successfully adding an IDTNode to the IDT
    */
   void increaseGlobalIDTNodeIndex()  { _nextIdx ++; }

   /**
    * @brief Get the IDTNode using index.
    * Before using this method for accessing IDTNode, flattenIDT() must be called.
    *
    * @return the IDT node
    */
   TR::IDTNode *getNodeByGlobalIndex(int32_t index);

   /**
    * @brief Flatten all the IDTNodes into a list.
    */
   void flattenIDT();

   void print();

   private:
   TR::Compilation* comp() { return _comp; }

   TR::Compilation *_comp;
   TR::Region&  _region;
   int32_t _nextIdx;
   uint32_t _totalCost;
   TR::IDTNode* _root;
   TR::IDTNode** _indices;
   };

/**
 * Accessing IDTNode by the priority of its cost and benefit.
 * This priority queue is to prepare an ordered list of inlining options as described in following page:
 * https://patents.google.com/patent/US10055210B2/en
 *
 * It has 3 private variables:
 * _idt: same as the input IDT.
 * _entries: a vector used to store the result IDTNodes which were popped out from _pQueue.
 * _pQueue: a priority queue keeps IDTNode in the order of priority from highest to lowest.
 *
 * The queue has overwrite the operator for comparing he IDTNodes:
 * - the node with the lowest benefit will have the highest priority;
 * - the queue breaks the tie by comparing the cost (node with lowest cost will have higher priority).
 *
 * Right after the priority queue is constructed, only the root of given IDT will be pushed into _pQueue.
 * When getting the IDTNode at the specific index at the priority queue, it will check if the node has
 * already been saved in _entries first. If yes, return the IDTNode at the given index from _entries, else:
 *    1. pop the IDTNode which has the highest priority from _pQueue;
 *    2. append the popped node at the end of _entries;
 *    3. push all children of the popped node into _pQueue;
 *    4. return the IDTNode at the given index from _entries.
 *
 * This priority queue traverses all predecessors in preorder before accessing a successor.
 */
class IDTPriorityQueue
   {
   public:
   IDTPriorityQueue(TR::IDT* idt, TR::Region& region);
   uint32_t size() { return _idt->getNumNodes(); }

   TR::IDTNode* get(uint32_t index);

   private:
   struct IDTNodeCompare
      {
      bool operator()(TR::IDTNode *left, TR::IDTNode *right)
         {
         TR_ASSERT_FATAL(left && right, "Comparing against null");
         if (left->getBenefit() == right->getBenefit()) return left->getCost() >= right->getCost();
         else return left->getBenefit() > right->getBenefit();
         }
      };

   typedef TR::vector<IDTNode*, TR::Region&> IDTNodeVector;
   typedef std::priority_queue<IDTNode*, IDTNodeVector, IDTNodeCompare> IDTNodePriorityQueue;

   TR::IDT* _idt;
   IDTNodePriorityQueue _pQueue;
   IDTNodeVector _entries;
   };
}

#endif
