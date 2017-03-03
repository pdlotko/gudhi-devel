/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef DOC_GUDHI_STAT_H_
#define DOC_GUDHI_STAT_H_

namespace Gudhi {

namespace Gudhi_stat {

/**  \defgroup Gudhi_stat Gudhi_stats
 *
 * \author   Pawel Dlotko
 *
 * @{
 *\section{Idea}

 *In order to perform most of the statistical tests and machine learning algorithms on a data one need to be able to perform only a very limited number of operations on them. Let us fix a representation of data of a type $\mathcal{A}$. To perform most of the statistical and machine learning operations one need to be able to compute average of objects of type $\mathcal{A}$ (so that the averaged object is also of a type $\mathcal{A}$), to compute distance between objects of a type $\mathcal{A}$, to vectorize object of a type $\mathcal{A}$ and to compute scalar product of a pair objects of a type $\mathcal{A}$.
 *
 *To put this statement into a context, let us assume we have two collections $c_1,...,c_n$ and $d_1,...,d_n$ of objects of a type $\mathcal{A}$. We want to verify if the average of those two collections are different by performing a permutation test.
 *First of all, we compute averages of those two collections: $C =$ average of $c_1,...,c_n$ and $D =$ average of $d_1,...,d_n$. Note that both $C$ and $D$ are of a type $\mathcal{A}$. Then we compute $d(C,D)$, a distance between $C$ and $D$.
 *Later we put the two collections into one bin:\\
 *\[B = { c_1,...,c_n,d_1,...,d_n }\]
 *Then we shuffle $B$, and we divide the shuffled version of $B$ into two classes: $B_1$ and $B_2$ (in this case, of the same cardinality). Then we compute averages $\hat{B_A}$ and $\hat{B_2}$ of elements in $B_1$ and $B_2$. Note that again, $\hat{B_1}$ and $\hat{B_2}$ are of a type $\mathcal{A}$. Then we compute their distance $d(\hat{B_1},\hat{B_2})$. The procedure of shuffling and dividing the set $B$ is repeated $N$ times (where $N$ is reasonably large number).
 *Then the p-value of a statement that the averages of  $c_1,...,c_n$ and $d_1,...,d_n$ is approximated by the number of times $d(\hat{B_1},\hat{B_2}) > d(C,D)$ divided by $N$.
 *
 *The permutation test reminded above can be performed for any type $\mathcal{A}$ which can be averaged, and which allows for computations of distances. 
 *
 *The Gudhi\_stat contains a collection of various representations of persistent homology that implements various concepts described below:
 *\enumerate{
 *\item{1}{Interface of a representation of persistence that allows averaging (so that the average object is of the same type).}
 *\item{2}{Interface of representation of persistence that allows computations of distances.}
 *\item{3}{Interface of representation of persistence that allows computations of scalar products.}
 *\item{4}{Interface of representation of persistence that allows vectorization.}
 *\item{5}{Interface of representation of persistence that allows computations of real--valued characteristics of objects.}
 *\}
 *
 *At the moment an implementation of the following representations of persistence are available (further details of those representations will be discussed later):
 *\enumerate{
 *\item{1}{Exact persistence landscapes (allow averaging, computation of distances, scalar products, vectorizations and real value characteristics).}
 *\item{2}{Persistence landscapes on a grid (allow averaging, computation of distances scalar products, vectorizations and real value characteristics).}
 *\item{3}{Persistence heat maps – various representations where one put some weighted or not Gaussian kernel for each point of diagram (allow averaging, computation of distances, scalar products, vectorizations and real value characteristics).}
 *\item{4}{Persistence vectors (allow averaging, computation of distances, scalar products, vectorizations and real value characteristics).}
 *\item{5}{Persistence diagrams / barcodes (allow computation of distances, vectorizations and real value characteristics).}
 *}
 *
 *Note that at the while functionalities like averaging, distances and scalar products are fixed, there is no canonical way of vectorizing and computing real valued characteristics of objects. Therefore the vectorizations and computation of real value characteristics procedures are quite likely to evolve in the furthering versions of the library. 
 *
 *The main aim of this implementation is to be able to implement various statistical methods, both on the level of C++ and on the level of python. The methods will operate on the functionalities offered by concepts. That means that the statistical and ML methods will be able to operate on any representation that implement the required concept (including the ones that are not in the library at the moment). That gives provides a framework, that is very easy to extend, for topological statistics.
 *
 *Below we are discussing the representations which are currently implemented in Gudhi\_stat:
 *
 *\section{Persistence Landscapes}
 *\label{sec:persistence_landscapes}
 *Persistence landscapes were originally proposed by Bubenik in~\cite{landscapes1}. Efficient algorithms to compute them rigorously were proposed by Bubenik and D{\l}otko in~\cite{landscapes2}. The idea of persistence landscapes is shortly summarized in below.
 *
 *To begin with, suppose we are given a point~$(b,d) \in \mathbb{R}^2$ in a
 *persistence diagram. With this point, we associate a piecewise
 *linear function~$f_{(b,d)} : \mathbb{R} \rightarrow [0,\infty)$, which is
 *defined as
 *
 *\begin{equation} \label{eq:basicLand}
 *  f_{(b,d)}(x) =
 *  \left\{ \begin{array}{ccl}
 *            0     & \mbox{ if } & x \not\in (b, d) \; , \\[2ex]
 *            x - b & \mbox{ if } & x \in \left( b, \frac{b+d}{2}
 *              \right] \; , \\[2ex]
 *            d - x & \mbox{ if } & x \in \left(\frac{b+d}{2},
 *              d \right) \; .
 *  \end{array} \right.
 *\end{equation}
 *
 *A \emph{persistence landscape} of the birth-death
 *pairs~$(b_i , d_i)$, where $i = 1,\ldots,m$, which constitute the given
 *persistence diagram is the sequence of functions $\lambda_k : \mathbb{R}
 *\rightarrow [0,\infty)$ for $k \in \mathbb{N}$, where~$\lambda_k(x)$
 *denotes the $k^{\rm th}$ largest value of the numbers~$f_{(b_i,d_i)}(x)$,
 *for $i = 1, \ldots, m$, and we define $\lambda_k(x) = 0$ if $k > m$.
 *Equivalently, this sequence of functions can be combined into a single
 *function $L : \mathbb{N} \times \mathbb{R} \to [0,\infty)$ of two
 *variables, if we define $L(k,t) = \lambda_k(t)$.
 *
 *The detailed description of algorithms used to compute persistence landscapes can be found in~\cite{landscapes2}. Note that this implementation provides exact representation of landscapes. That have many advantages, but also a few drawbacks. For instance, as discussed in~\cite{landscapes2}, the exact representation of landscape may be of quadratic size with respect to the input persistence diagram. It may therefore happen that, for very large diagrams, using this representation may be memory--prohibitive. In such a case, there are two possible ways to proceed:
 *\enumerate{
 *\item{1}{Use non exact representation on a grid described in the Section~\ref{sec:landscapes_on_grid}}.
 *\item{2}{Compute just a number of initial nonzero landscapes. This option is available from C++ level}. 
 *}
 * 
 * 
 *\section{Persistence Landscapes on a grid}
 *\label{sec:landscapes_on_grid}
 *This is an alternative, not--exact, representation of persistence landscapes defined in the Section~\ref{sec:persistence_landscapes}. Unlike in the Section~\ref{sec:persistence_landscapes} we build a representation of persistence landscape by sampling its values on a finite, equally distributed grid of points. Since, the persistence landscapes that originate from persistence diagrams have slope $1$ or $-1$, we have an estimate of a region between the grid points where the landscape cab be located. That allows to estimate an error make when performing various operations on landscape. Note that for average landscapes the slope is in range $[-1,1]$ and similar estimate can be used. 
 *
 *Due to a lack of rigorous description of the algorithms to deal with this non--rigorous representaion of persistence landscapes in the literature, we are providing a short discussion of them in below.
 *
 *Let us assume that we want to compute persistence landscape on a interval $[x,y]$. Let us assume that we want to use $N$ grid points for that purpose. Then we will sample the persistence landscape on points $x_1 = x , x_2 = x + \frac{y-x}{N}, \ldots , x_{N} = y$. Persistence landscapes are represented as a vector of vectors of real numbers. Assume that i-th vector consist of $n_i$ numbers sorted from larger to smaller. They represent the values of the functions $\lambda_1,\ldots,\lambda_{n_i}$ ($\lambda_{n_i+1}$ and the functions with larger indices are then zero functions) on the i-th point of a grid, i.e. $x + i \frac{y-x}{N}$. 
 *
 *When averaging two persistence landscapes represented by a grid we need to make sure that they are defined in a compatible grids. I.e. the intervals $[x,y]$ on which they are defined are the same, and the numbers of grid points $N$ are the same in both cases. If this is the case, we simply compute point-wise averages of the entries of corresponding vectors\footnote{In this whole section we assume that if one vector of numbers is shorter than another, we extend the shorter one with zeros so that they have the same length.}
 *
 *Computations of distances between two persistence landscapes on a grid is not much different than in the rigorous case. In this case, we sum up the distances between the same levels of corresponding landscapes. For fixed level, we approximate the landscapes between the corresponding constitutive points of landscapes by linear functions, and compute the $L^p$ distance between them.
 *
 *Similarly as in case of distance, when computing the scalar product of two persistence landscapes on a grid, we sum up the scalar products of corresponding levels of landscapes. For each level, we assume that the persistence landscape on a grid between two grid points is approximated by linear function. Therefore to compute scalar product of two corresponding levels of landscapes, we sum up the integrals of products of line segments for every pair of constitutive grid points. 
 *
 *Note that for this representation we need to specify a few parameters:
 *\enumerate{
 *\item{1}{Begin and end point of a grid -- the interval $[x,y]$ (real numbers).}
 *\item{2}{Number of points in a grid (positive integer $N$).}
 *}
 *
 *Note that the same representation is used in TDA R-package~\cite{tda}.
 *
 *\section{Persistence heat maps}
 *This is a general class of discrete structures which are based on idea of placing a kernel in the points of persistence diagrams. This idea appeared in work by many authors over the last 15 years. As far as we know this idea was firstly described in the work of Bologna group in~\cite{bologna1} and~\cite{bologna2}. Later it has been described by Colorado State University group in~\cite{Henry}. The presented paper in the first time provide a discussion of stability of the representation. Also, the same ideas are used in construction of two recent kernels used for machine learning:~\cite{yasu} and~\cite{uli}. Both the kernel's construction uses interesting ideas to ensure stability of the representation with respect to Wasserstein metric. In the kernel presented in~\cite{yasu}, a scaling function is used to multiply the Gaussian kernel in the way that the points close to diagonal got low weight and consequently do not have a big influence on the resulting distribution. In~\cite{uli} for every point $(b,d)$ two Gaussian kernels are added: first, with a weight $1$ in a point $(b,d)$, and the second, with the weight $-1$ for a point $(b,d)$. In both cases, the representations are stable with respect to 1-Wasserstein distance.
 *
 *In Gudhi\_stat we currently implement a discretization of the distributions described above. The base of this implementation is 2-dimensional array of pixels. Each pixel have assigned a real value which is a sum of values of distributions induced by each point of the persistence diagram. At the moment we compute the sum of values on a center of a pixels. It can be easily extended to any other function (like for instance sum of integrals of the intermediate distribution on a pixel). 
 *
 *The parameters that determine the structure are the following:
 *\enumerate{
 *\item{1}{A positive integer $k$ determining the size of the kernel we used (we always assume that the kernels are square).}
 *\item{2}{A filter: in practice a square matrix of a size $2k+1 \times 2k+1$. By default, this is a discretization of N(0,1) kernel.}
 *\item{3}{The box $[x_0,x_1]\times [y_0,y_1]$ bounding the domain of the persistence image.}
 *\item{4}{Scaling function. Each Gaussian kernel at point $(p,q)$ gets multiplied by the value of this function at the point $(p,q)$.}
 *\item{5}{A boolean value determining if the space below diagonal should be erased or not. To be precise: when points close to diagonal are given then sometimes the kernel have support that reaches the region below the diagonal. If the value of this parameter is \emph{true}, then the values below diagonal can be erased.}
 *}
 *
 *\section{Persistence vectors}
 *This is a representation of persistent homology in a form of a vector which was designed for an application in 3d graphic in~\cite{vectors}. Below we provide a short description of this representation.
 *
 *Given a persistence diagram $D = \{ (b_i,d_i) \}$, for every pair of birth--death points $(b_1,d_1)$ and $(b_2,d_2)$ we compute the following three distances:
 *\enumerate{
 *\item{1}{$d( (b_1,d_1) , (b_2,d_2) )$.}
 *\item{2}{$d( (b_1,d_1) , (\frac{b_1,d_1}{2},\frac{b_1,d_1}{2}) )$.}
 *\item{3}{$d( (b_2,d_2) , (\frac{b_2,d_2}{2},\frac{b_2,d_2}{2}) )$.}
 *\}
 *We pick the smallest of those and add it to a vector. The obtained vector of numbers is then sorted in decreasing order. This way we obtain a \emph{persistence vector} representing the diagram.
 *
 *Given two persistence vectors, the computation of distances, averages and scalar products is straightforward. Average is simply a coordinate-wise average of a collection of vectors. In this section we assume that the vectors are extended by zeros if they are of a different size. To compute distances we compute absolute value of differences between coordinates. A scalar product is a sum of products of values at the corresponding positions of two vectors. 
 *
 * \copyright GNU General Public License v3.
 */
/** @} */  // end defgroup Gudhi_stat

}  // namespace cubical_complex

namespace Gudhi_stat = Gudhi_stat;

}  // namespace Gudhi

#endif  // DOC_GUDHI_STAT_H_
