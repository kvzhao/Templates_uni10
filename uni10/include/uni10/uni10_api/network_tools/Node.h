/****************************************************************************
*  @file uni10/uni10_api/node.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Generic header file for uni10 data strucutres
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_NODE_H__
#define __UNI10_NODE_H__

#include "uni10/uni10_api/Bond.h"

namespace uni10{

  typedef struct{
    int b1;
    int b2;
  }_Swap;

  class Bond;

  template<typename uni10_type>
    class UniTensor;

  template <typename uni10_type>
    class Node {

      public:
        Node();
        Node(UniTensor<uni10_type>* Tp);
        Node(const Node& nd);
        Node(std::vector<Bond>& _bonds, std::vector<uni10_int32>& _labels);
        ~Node();
        Node contract(Node* nd);
        uni10_float32 metric(Node* nd);

        friend std::ostream& operator<< (std::ostream& os, const Node& nd);

        template<typename _uni10_type>
          friend class Network;

      private:
        UniTensor<uni10_type>* T;   //if T != NULL, it is leaf node
        std::vector<uni10_int32> labels;
        std::vector<Bond> bonds;
        uni10_uint64 elemNum;
        std::string name;
        Node* parent;
        Node* left;
        Node* right;
        uni10_float32 point;
        uni10_uint64 cal_elemNum(std::vector<Bond>& _bonds);
        void delink();

    };

}; /* namespace uni10 */

#endif
