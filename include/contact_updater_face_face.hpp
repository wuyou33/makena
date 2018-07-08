#ifndef _MAKENA_CONTACT_UPDATER_FACE_FACE_HPP_
#define _MAKENA_CONTACT_UPDATER_FACE_FACE_HPP_

#include <memory>
#include <array>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <exception>
#include <stdexcept>
#include <cmath>
#include <cstdarg>
#include <iomanip> 

#include "convex_rigid_body.hpp"
#include "binary_dilation.hpp"
#include "contact_pair_info.hpp"
#include "loggable.hpp"
#include "intersection_convex_polygon_2d.hpp"

#ifdef UNIT_TESTS
#include "gtest/gtest_prod.h"
#endif


/**
 * @file contact_updater_face_face.hpp
 *
 * @brief assuming the current contact feature pair are face-face,
 *        check if it needs to be updated to another feature pair
 *        based on the current face normal directions and intersection of
 *        the two faces.
 *       
 *  Legend:
 *
 *          //{eps} : Parallel within angular tolerance eps
 *
 *         _|_{eps} : Perpendicular within angular tolerance eps
 *
 *          Nx : Normal vector. Usually a face normal.
 *
 *          Dx : Direction vector. Usually half edge direction from src to dir.
 *
 *          ·  : dot product
 *
 *          x  : cross product
 *
 *         |V| : Vector norm 2
 *
 *
 *  FACE1-FACE2
 *      
 *      Let Nf1 be the normal of FACE1
 *      Let Nf2 be the normal of FACE2
 *
 *      - If Nf1 //{2.0eps}  => FACE1-FACE2 (NO CHANGE)
 *
 *      - Else
 *
 *        Let D12 be the unit tilting direction vector norm(Nf1 + Nf2)
 *
 *        Let INT be the 2D convex intersection polygon of FACE1 and FACE 2
 *        on the XY-plane perpendicular to norm(Nf1 - Nf2).
 *        Find the extremal vertex Pint or edge (Pint1, Pint2) 
 *        of INT along D12. 
 *        If it's an edge (Pint1, Pint2), then (Pint1,Pint2) _|_{eps} D12.
 *        If it's a vertex Pint, then 
 *
 *            NOT [(PintPrev, Pint    ) _|_{eps} D12] &&
 *            NOT [(Pint,     PintNext) _|_{eps} D12] 
 *
 *        - If it is a vertex Pint:
 *
 *          - If Pint = IT_EDGE_EDGE (Eint1,Eint2) => Eint1-Eint2 (EDGE-EDGE)
 *
 *                               \
 *                                *                ^
 *                                |Eint1          /
 *                                |              /
 *                        o-------+-----o      D12
 *                       /.Eint2..|      \
 *                       .........|
 *                      ..........*
 *                      ........./
 *
 *          - If Pint = IT_VERTEX_INTERIOR (Vint1, FACE2)
 *                                             => Vint1-FACE2 (VERTEX-FACE)
 *
 *                              Vint1            ^
 *                       *------*               /
 *                      Eint1...|              /
 *                       .......|           D12
 *                       .......|Eint2      
 *                              *
 *
 *          - If Pint = IT_INTERIOR_VERTEX (FACE1, Vint2)
 *                                             => FACE1-Vint2 (FACE-VERTEX)
 *
 *                              Vint2            ^
 *                       *------*               /
 *                      Eint1...|              /
 *                       .......|           D12
 *                       .......|Eint2      
 *                              *
 *
 *
 *          - If Pint = IT_EDGE_VERTEX(Eint1, Vint2)
 *            There are many possible configurations of edges and vertices
 *            as follows. However, they all end up in the same update, 
 *            which is EDGE-VERTEX.
 *
 *              o  *       o  *        *  o       *o          *o  
 *               \ |      ..\ |    ....| /        ||TP     ...||
 *                \|      ...\|    ....|/         ||       ...||
 *                 oTP    ....o    ....o          |o       ...|o
 *                /|      .../|    .../|         /|        ../|
 *               / |      ../ |    ../ |        / |        ./ |
 *              o  *       o  *     o  *       o  *        o  *
 *
 *
 *              o  *       o  *      o  *        *o         *o
 *              .\ |        \ |      .\ |        ||TP    ...||
 *              ..\|         \|      ..\|        ||      ...||
 *              ...o          o      ...o        |o      ...|o
 *              ...|\         ||TP   ...|        ||TP    ...||
 *              ...| \        ||     ...|        ||      ...||
 *                 *  o       *o        *o       *o         *o
 *
 *                                    ==> Eint1-Vint2  (EDGE-VERTEX)
 *
 *         - If Pint = IT_VERTEX_EDGE(Vint1, Eint2)
 *                                    ==> Vint1-Eint2  (VERTEX-EDGE)
 *
 *
 *
 *         - If Pint = VERTEX-VERTEX (Vint1, Vint2)
 *            There are many possible configurations of edges and vertices
 *            as follows. However, they all end up in the same update, 
 *            which is VERTEX-VERTEX.
 *                                    ==> Vint1-Vint2  (VERTEX-VERTEX)
 *
 *              *          *          *         *          *o
 *               \          \          \         \         ..\
 *            *-_ \TP    o-_ \      o-_ \         \TP      ...\
 *              _-*o     .._-*o     .._-*o    *o==*o       *--*o
 *            o-  /      *-  /      o-  /         /           /
 *               /          /          /         /           /
 *              o          o          *         o           o
 *
 *              *o         *          *        *     *o
 *             ..\          \          \        \    .\\
 *             ...\          \          \        \   ..\\
 *             o--*o     *o==*oTP    o--*o    o--*o  ...*o
 *                /          /       .../     .../   ..//
 *               /          /        ../      ../    .//
 *              *          o          o*       *o    *o
 *
 *
 *        - If it is an edge (Pint1, Pint2)
 *
 *          - If Pint1 = IT_EDGE_EDGE (Eint11, Eint12)
 *            - If Pint2 = IT_EDGE_EDGE (Eint21, Eint22)
 *
 *              - Eint11==Eint21 on FACE1 => Eint11-FACE2 (EDGE-FACE)
 *
 *                D12 -->                      D12 -->
 *
 *                   o  *                       *  o
 *                ....\ |Eint11==Eint21         | /Eint12
 *               Eint12\|                       |/
 *                ......+                       +
 *                ......|\                     /|
 *                ......| \                   /.|
 *                ......|  o                 o..|
 *                ......|  |                 |..|Eint11==Eint21
 *                ......|  |                 |..|
 *                ......|  o                 o..|
 *                ......| /                   \.|
 *                ......|/                     \|
 *                ......+                       +
 *               Eint22/|                       |\
 *                ..../ |                       | \Eint22
 *                   o  *                       *  o
 *
 *              - Eint12==Eint22 on FACE2 => FACE1-Eint12 (FACE-EDGE)
 *
 *                D12 -->                      D12 -->
 *
 *                   *  o                       o  *
 *                ....\ |Eint12==Eint22         | /Eint11
 *               Eint11\|                       |/
 *                ......+                       +
 *                ......|\                     /|
 *                ......| \                   /.|
 *                ......|  *                 *..|
 *                ......|  |                 |..|Eint12==Eint22
 *                ......|  |                 |..|
 *                ......|  *                 *..|
 *                ......| /                   \.|
 *                ......|/                     \|
 *                ......+                       +
 *               Eint21/|                       |\
 *                ..../ |                       | \Eint21
 *                   *  o                       o  *
 *
 *
 *            - If Pint2 = IT_EDGE_VERTEX (Eint21, Vint22)
 *
 *              - Vint22 incident to Eint12 => FACE1-Eint12 (FACE-EDGE)
 *
 *                D12 -->                      D12 -->
 *
 *                   .* o                      o *
 *                   ..\|                      |/Eint11
 *                   ...+                      +
 *                   ...|\Eint11              /|
 *                Eint12| *                  *.|
 *                   ...| |                  |.|Eint12
 *                   ...| *                  *.|
 *                   ...|/                    \|
 *                   o--oVint22          Vint22o--o
 *                     /Eint21                  \Eint21
 *                    *                          *
 *
 *              - Eint11==Eint21 on FACE1 => Eint11-FACE2 (EDGE-FACE)
 *
 *                D12 -->                      D12 -->
 *
 *                   o *                       * o
 *                  ..\|                       |/Eint12
 *                  ...+                       +
 *                  ...|\Eint12               /|
 *               Eint11| \                   /.|Eint11==Eint21
 *             ==Eint21|  o                 o..|
 *                  ...|  |                 |..|
 *                  ...|  o                 o..|
 *                  ...| /                   \.|
 *                  ...|/                     \|
 *                  o--oVint22           Vint22o--o
 *                     |                       |
 *                     *                       *
 *
 *            - If Pint2 = IT_VERTEX_EDGE (Vint21, Eint22)
 *
 *              - Vint21 incident to Eint11 => Eint11-FACE2 (EDGE-FACE)
 *
 *                D12 -->                      D12 -->
 *
 *                   .o *                      * o
 *                   ..\|                      |/Eint12
 *                   ...+                      +
 *                   ...|\Eint12              /|
 *                Eint11| o                  o.|
 *                   ...| |                  |.|Eint11
 *                   ...| o                  o.|
 *                   ...|/                    \|
 *                   *--*Vint21          Vint21*--*
 *                     /Eint22                  \Eint22
 *                    o                          *
 *
 *              - Eint12==Eint22 on FACE2 => FACE1-Eint12 (FACE-EDGE)
 *
 *                D12 -->                      D12 -->
 *
 *                   * o                       o *
 *                  ..\|                       |/Eint11
 *                  ...+                       +
 *                  ...|\Eint11               /|
 *               Eint12| \                   /.|Eint12==Eint22
 *             ==Eint22|  *                 *..|
 *                  ...|  |                 |..|
 *                  ...|  *                 *..|
 *                  ...| /                   \.|
 *                  ...|/                     \|
 *                  *--*Vint21           Vint21*--*
 *                     |                       |
 *                     o                       o
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *                                          => Eint11-FACE2 (EDGE-FACE)
 *              (Invariant: Vint21 is incident to Eint11.)
 *
 *                  D12 -->
 *
 *                 |..\
 *                 |...*Vint21
 *                 o...|
 *                  \..|
 *                   \.|Eint11
 *                    \|
 *                     +
 *                     |\
 *                     | \Eint12
 *                     *  o
 *
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *                                          => FACE1-Eint12 (FACE-EDGE)
 *              (Invariant: Vint22 is incident to Eint12.)
 *
 *                  D12 -->
 *
 *                 |..\
 *                 |...oVint22
 *                 *...|
 *                  \..|
 *                   \.|Eint12
 *                    \|
 *                     +
 *                     |\
 *                     | \Eint11
 *                     o  *
 *
 *
 *            - If Pint2 = IT_VERTEX__VERTEX (Vint21, Vint22)
 *
 *              - Vint21 is incident to Eint11 => Eint11-FACE2 (EDGE-FACE)
 *
 *                  D12 -->               D12 -->
 *
 *                   o  *                     *  o
 *              Eint12\ |Eint11         Eint11| /Eint12
 *                   ..\|                     |/
 *                   ...+                     +
 *                   ...|\                   /|
 *                   ...| \                 /.|
 *                   ...|  o               o..|
 *                   ...|  |               |..|
 *                   ...|  |               |..|
 *                   ...|  o               o..|
 *                   ...| /                 \.|
 *                   ...|/                   \|
 *                   *--*oVint21==     Vint21o*--*
 *                     /  Vint22     ==Vint22  \
 *                    o                         o
 *
 *
 *              - Vint22 is incident to Eint12 => FACE1-Eint12 (FACE-EDGE)
 *
 *                  D12 -->               D12 -->
 *
 *                   *  o                     o  *
 *              Eint11\ |Eint12         Eint12| /Eint11
 *                   ..\|                     |/
 *                   ...+                     +
 *                   ...|\                   /|
 *                   ...| \                 /.|
 *                   ...|  *               *..|
 *                   ...|  |               |..|
 *                   ...|  |               |..|
 *                   ...|  *               *..|
 *                   ...| /                 \.|
 *                   ...|/                   \|
 *                   o--*oVint21==     Vint21*o--o
 *                     /  Vint22      ==Vint22 \
 *                    *                         *
 *
 *
 *          - If Pint1 = IT_EDGE_VERTEX (Eint11, Vint12)
 *
 *            - If Pint2 = IT_EDGE_VERTEX (Eint21, Vint22) 
 *
 *              - If (Vint12, Vint22) is on FACE2, 
 *
 *                - If Eint11 == Eint21  => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *                                  
 *                    D12 -->               D12 -->
 *
 *                       \                   /
 *                      o *               o *
 *                     ..\|                \|
 *                     ...o                 o
 *                     ...|                 |
 *                     ...o                 o
 *                     ../|                /|
 *                      o *               o *
 *                       /                   \
 *
 *                - If NOT ((Vint12,Vint22) //{eps} Eint11 ) &&
 *                     NOT ((Vint12,Vint22) //{eps} Eint21 )
 *
 *                                    => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *
 *                    D12 -->               D12 <--
 *          
 *                     *  o                 *     o
 *                     .\ |                  \   /
 *                     ..\|                   \ /
 *                     ...oVint12              oVint12
 *                     ...|\Eint11             |\Eint11
 *                     ...| *                  |.*
 *                     ...| |                  |.|
 *                     ...| *                  |.*
 *                     ...|/Eint21             |/Eint21
 *                     ...oVint22              oVint22
 *                     ../|                   / \
 *                     ./ |                  /   \
 *                     *  o                 *     o
 *
 *                - If ((Vint12,Vint22) //{eps} Eint11 )
 *
 *                                    => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                    D12 -->               D12 -->
 *
 *                        \                 \
 *                    ..o  *                 *  o
 *                    ...\ |Eint11     Eint11| /
 *                  Vint12o|                 |oVint12
 *                    ....||                 ||
 *                    ....|*                 *|
 *                    ....||                 ||
 *                  Vint22o|                 |oVint22
 *                    .../ |Eint21     Eint21| \
 *                    ..o  *                 *  o
 *                        /                 /
 *
 *                - If ((Vint12,Vint22) //{eps} Eint12 )
 *
 *                                    => Eint12-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                    D12 -->               D12 -->
 *
 *                        \                 \
 *                    ..o  *                 *  o
 *                    ...\ |Eint11     Eint11| /
 *                  Vint12o|                 |oVint12
 *                    ....||                 ||
 *                    ....|*                 *|
 *                    ....||                 ||
 *                  Vint22o|                 |oVint22
 *                    .../ |Eint21     Eint21| \
 *                    ..o  *                 *  o
 *                        /                 /
 *
 *              - If Eint11 == Eint21 && (Vint12, Vint22) is not on FACE2
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                    D12 -->                     D12 --> 
 *                      \                           \
 *                    o  *                           *  o
 *                    .\ |Eint11==Eint21             | /
 *                    ..\|                           |/
 *                    ...oVint12                     o
 *                    ...|\                         /|
 *                    ...| o                       o.|
 *                    ...| |                       |.|
 *                    ...| o                       o.|
 *                    ...|/                         \|
 *                    ...oVint22                     o
 *                    ../|                           |\
 *                    ./ |                           | \
 *                    o  *                           *  o
 *                      /                           /
 *
 *
 *            - If Pint2 = IT_VERTEX_EDGE (Vint21, Eint22)
 *
 *              - If Vint21 is incident to Eint11 &&
 *                   Vint12 is incident to Eint22
 *                (Invariant:  Eint11 //{eps} Eint22)
 *
 *                                               => Eint11-Eint22 (EDGE-EDGE)
 *
 *                  D12 -->             D12 -->
 *
 *                   \                   
 *                    *                   *  o
 *                ..\ |Eint11       Eint11| /
 *                ...\|                   |/
 *                ....oVint12             oVint12
 *                ....|                   ||
 *                ....|                   ||
 *              Eint11|                   ||
 *                ....|                   ||
 *                ....|                   ||
 *                ... *Vint21             *Vint21
 *                .../|                  /|
 *                ../ |Eint22           / |Eint22
 *                 *  o                *  o
 *
 *
 *              - If Vint21 is incident to Eint11
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                  D12 -->             D12 -->
 *
 *                    *                   *
 *                    |                   |
 *                 o--oVint12       Vint12o--o
 *                ....|\                 /|
 *                ....| \               /.|
 *              Eint11|  o             o..|Eint11
 *                ....| /Eint22   Eint22\.|
 *                ....|/                 \|
 *                ... *Vint21       Vint21*
 *                .../|                   |\ 
 *                  o |                   | o
 *                    *                   *
 *
 *
 *              - If Vint12 is incident to Eint22
 *
 *                                               => FACE1-Eint22 (FACE-EDGE)
 *
 *                  D12 -->             D12 -->
 *
 *                    o                   o
 *                    |                   |
 *                 *--*Vint21       Vint21*--*
 *                ....|\                 /|
 *                ....| \               /.|
 *              Eint22|  *             *..|Eint22
 *                ....| /Eint11   Eint11\.|
 *                ....|/                 \|
 *                ... oVint12       Vint12o
 *                .../|                   |\ 
 *                  * |                   | *
 *                    o                   o
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *                  D12 -->
 *
 *                  o  *
 *                 ..\ |Eint11
 *                 ...\|
 *                 ....oVint12
 *                 ....|\
 *                 ....| \
 *               Vint21*  \
 *                 .../    o
 *                   *
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *              (Invariant: Edge (Vint12,Vint22) is on FACE2)
 *                                        => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *                   D12 -->
 *
 *                   o     *             
 *                    \   /Eint11
 *                     \ /
 *                      oVint12
 *                     /||
 *                    /.||
 *                   /..||
 *                  /...oVint22
 *                 *  ./
 *                    /
 *                   o
 *
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *              - If Edge (Vint12,Vint22) is on FACE2 &&
 *                        Vint21 is incident to Eint11
 *                                      => Eint11-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                 D12 -->                  D12 -->
 *
 *                     *  o                 o  *
 *               Eint11| /                  .\ |Eint11
 *                     |/                   ..\|
 *                     |oVint12             ..o|Vint12
 *                     ||                   ..||
 *                     ||                   ..||
 *                     ||                   ..||
 *                     *oVint21==Vint22     ..o*Vint21==Vint22
 *                    /  \                  ./|
 *                   *    o                 * o
 *
 *              - If Edge (Vint12,Vint22) is on FACE2
 *                                        => FACE1-(Vint12,Vint22) (FACE-EDGE)
 *
 *                   D12 <--               D12 -->
 *
 *                  *     o                  *
 *             Eint11\   /                    \Eint11
 *                    \ /                      \
 *               Vint12o                    o---oVint12
 *                    ||\                   ....|\
 *                    ||.\                  ....| \
 *                    ||..*                 ....|  *
 *                    ||..|                 ....|  |
 *                    ||..|                 ....|  |
 *                    ||..*                 ....|  *
 *                    ||./                  ....| /
 *                    ||/                   ....|/
 *                     o*Vint21==Vint22     ...o*Vint21==Vint22
 *                    / \                   ../|
 *                   /   \                  ./ |
 *                  *     o                 o  *
 *
 *              - If Vint21 is incident to Eint11
 *                                               => Eint11-FACE2 (EDGE-FACE)
 *
 *                   D12 -->            D12 -->
 *
 *                  o  *                       *  o
 *                  .\ |Eint11                 | /
 *                  ..\|                       |/
 *                  ...oVint12                 oVint12
 *                  ...|\                     /|
 *                  ...| \                   /.|
 *                  ...|  o                 o..|
 *                  ...|  |                 |..|Eint11
 *                  ...|  |                 |..|
 *                  ...|  o                 o..|
 *                  ...| /                   \.|
 *                  ...|/                     \|
 *                  o--*oVint21==Vint22        *oVint21==Vint22
 *                    /                       /
 *                   /                       /
 *                  *                       *
 *
 *
 *          - If Pint1 = IT_VERTEX_EDGE (Vint11, Eint12)
 *
 *            - If Pint2 = IT_VERTEX_EDGE (Vint21, Eint22)
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *
 *          - If Pint1 = IT_VERTEX_INTERIOR (Vint11, *)
 *
 *            - If Pint2 = IT_VERTEX_INTERIOR (Vint21, *)
 *              (Invariant: Edge(Vint11,Vint21) is on FACE1)
 *                                     ==> (Vint11,Vint21)-FACE2  (EDGE-FACE)
 *               D12 -->
 *
 *                *
 *                .\
 *                ..\
 *                ...*Vint11
 *                ...|
 *                ...|
 *                ...|
 *                ...*Vint21
 *                ../
 *                ./
 *                *
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *
 *                FORBIDDEN
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *                (Invariant: Edge (Vint11,Vint21) is on FACE1)
 *                                     ==> (Vint11,Vint21)-FACE2  (EDGE-FACE)
 *
 *                  D12 -->
 *
 *                           o
 *                 ...\     /
 *                 ....*Vint11
 *                 ....|  /
 *                 ....| /
 *                 ....|/
 *                 o--*oVint21==Vint22
 *                    /
 *                   *
 *
 *          - If Pint1 = IT_INTERIOR_VERTEX (*, Vint12)
 *
 *            - If Pint2 = IT_INTERIOR_VERTEX (*, Vint22)
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *          - If Pint1 = IT_VERTEX_VERTEX (Vint11, Vint12)
 *
 *            - If Pint2 = IT_VERTEX_VERTEX (Vint21, Vint22)
 *
 *                           ==> (Vint11,Vint21)-(Vint12,Vint22) (EDGE-EDGE)
 *
 *                D12 -->     D12 -->         D12 -->
 *
 *                *    o         *              o
 *                 \  /           \              \
 *                  \/             \              \
 *                  *o          o--o*         *---*o
 *                  ||          ...||          ...||
 *                  ||          ...||          ...||
 *                  ||          ...||          ...||
 *                  *o          o--o*         o---o*
 *                  /\             /              /
 *                 /  \           /              /
 *                *    o         *              *
 */

namespace Makena {

using namespace std;

using IntSec2D = IntersectionFinderConvexPolygon2D;

class ContactUpdater_FACE_FACE : public Loggable {

  public:

    ContactUpdater_FACE_FACE(
        ConvexRigidBody&         body1,
        ConvexRigidBody&         body2,
        ContactPairInfo&         info,
        const double&            epsilonZero,
        const double&            epsilonAngle,
        std::ostream&            logStream
    );

    ~ContactUpdater_FACE_FACE();

    /** @brief main function
     *
     *  @return true if it needs to run GJK to find the contact pair.
     */
    bool update();

#ifdef UNIT_TESTS
  public:
#else
  private:
#endif

    void checkForTilt(bool& aligned, bool& abort);

    Mat3x3 rotMatAlignZDirToZAxis(const Vec3& zDir);

    void rotateWorld(const Vec3& zDir);

    bool findIntersection2D();

    Vec2 findTilt();

    void findExtremeFeature2D(const Vec2& dir, long& index1, long& index2);

    bool isVertexIncidentToEdge(VertexIt v, HalfEdgeIt h);

    bool findHalfEdgeBody1(
                  const long  index1, const long  index2, HalfEdgeIt& heOut);

    bool findHalfEdgeBody2(
                  const long  index1, const long  index2, HalfEdgeIt& heOut);

    bool areParallel(
        EdgeIt     eit1,
        const long body1,
        EdgeIt     eit2,
        const long body2
    );

    void dispatchAndUpdate(const long index1, const long index2);

    long edgeIndexA(IntSec2D::OutputElem& e);
    long edgeIndexB(IntSec2D::OutputElem& e);

    void processIntsec_EDGE_EDGE       (const long index);
    void processIntsec_EDGE_VERTEX     (const long index);
    void processIntsec_VERTEX_EDGE     (const long index);
    void processIntsec_VERTEX_VERTEX   (const long index);
    void processIntsec_VERTEX_INTERIOR (const long index);
    void processIntsec_INTERIOR_VERTEX (const long index);

    void processIntsec_EDGE_EDGE__EDGE_EDGE(
                                        const long index1, const long index2);
    void processIntsec_EDGE_EDGE__EDGE_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_EDGE_EDGE__VERTEX_EDGE(
                                        const long index1, const long index2);
    void processIntsec_EDGE_EDGE__VERTEX_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_EDGE_EDGE__VERTEX_INTERIOR(
                                        const long index1, const long index2);
    void processIntsec_EDGE_EDGE__INTERIOR_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_EDGE_VERTEX__EDGE_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_EDGE_VERTEX__VERTEX_EDGE(
                                        const long index1, const long index2);
    void processIntsec_EDGE_VERTEX__VERTEX_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_EDGE_VERTEX__VERTEX_INTERIOR(
                                        const long index1, const long index2);
    void processIntsec_EDGE_VERTEX__INTERIOR_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_EDGE__VERTEX_EDGE(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_EDGE__VERTEX_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_EDGE__VERTEX_INTERIOR(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_EDGE__INTERIOR_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_VERTEX__VERTEX_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_VERTEX__VERTEX_INTERIOR(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_VERTEX__INTERIOR_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_INTERIOR__VERTEX_INTERIOR(
                                        const long index1, const long index2);
    void processIntsec_VERTEX_INTERIOR__INTERIOR_VERTEX(
                                        const long index1, const long index2);
    void processIntsec_INTERIOR_VERTEX__INTERIOR_VERTEX(
                                        const long index1, const long index2);

    void checkAndUpdateEdgeCases();
    void updateNonTilting_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_EDGE();
    void updateNonTilting_EDGE_VERTEX();
    void updateNonTilting_EDGE_VERTEX_EDGE_VERTEX();
    void updateNonTilting_EDGE_VERTEX_VERTEX_EDGE();
    void updateNonTilting_EDGE_VERTEX_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_EDGE_EDGE_VERTEX();
    void updateNonTilting_VERTEX_EDGE_VERTEX_EDGE();
    void updateNonTilting_VERTEX_EDGE_VERTEX_VERTEX();
    void updateNonTilting_VERTEX_VERTEX_EDGE_VERTEX();
    void updateNonTilting_VERTEX_VERTEX_VERTEX_EDGE();
    void updateNonTilting_VERTEX_VERTEX_VERTEX_VERTEX();

    ConvexRigidBody&             mBody1;
    ConvexRigidBody&             mBody2;
    ContactPairInfo&             mInfo;
    const double                 mEpsilonZero;
    const double                 mEpsilonAngle;

    Mat3x3                       mRotMat1;
    Mat3x3                       mRotMat2;
    Vec3                         mCom1;
    Vec3                         mCom2;
    vector<VertexIt>             mVertices1;
    vector<VertexIt>             mVertices2;
    vector<HalfEdgeIt>           mHalfEdges1;
    vector<HalfEdgeIt>           mHalfEdges2;
    vector<IntSec2D::OutputElem> mIntsec;
    long                         mLastIndexA;
    long                         mLastIndexB;
};


}// namespace Makena

#endif /*_MAKENA_GJK_CONTACT_UPDATER_FACE_FACE_HPP_*/
