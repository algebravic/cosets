Maximum independent sets in $`Z^n/(2 D_n)`$
===========================================
Denote by $`D_n`$ the $n$-dimensiona sub-lattice of $\mathbb{Z}^n$:
$`\{x \in \mathbb{Z}^n: \sum_i x_n \equiv 0 \mod 2\}`$.  We are
interested in the quotient $`C_n := \mathbb{Z}^n/(2 D_n)`$.  Note that
since $`2 \mathbb{Z}^n \supset D_n`$, that $`C_n`$ is naturally a code
over the ring $\mathbb{Z}/4\mathbb{Z}$.  As such we have a natural
graph structure on $`C_n`$ induced by the Hamming metric.  The
question at hand is: what is the cardinality of the largest
independent set on $`C_n`$?

From Veit Elser:

For constructing the graph I like to work with a convenient choice of
coset representatives. Here’s the one I like:

1) First coordinate is an element of $`Z_4`$, say -1,0,1,2.
2) All other coordinates are elements of $`Z_2`$, say 0,1.

You get to this by first adding multiples of 4 to put all the
coordinates into $`Z_4`$. Next, you apply pairs of 2’s like (2,2,0,…)
or (2,0,2,…) or (2,-2,0,…) etc. so that all coordinates are 0 and 1
except that a single coordinate may also be -1 or 2. Finally, you
again add pairs of 2’s to move the -1 or 2 (if you have one) into the
first position.

Okay, that’s your $4 \cdot 2^{n-1} = 2^{n+1}$ points.

Here’s the rule for when to connect pairs of points.

Let point 1 have coordinates $`(z_1, b_1)`$, where $`z_1`$ is in
$`Z_4`$ and $`b_1`$ is a bit string, and similarly for point 2. Let
$`h(b_1,b_2)`$ be the Hamming distance between $`b_1`$ and $`b_2`$,
and $`g(z_1,z_2)`$ the $`Z_4`$ distance between $`z_1`$ and
$`z_2`$. For example, $g(-1,2)=1$. Points are connected if their
Euclidean distance is less than 2 (sum of squares less than 4). That
gives you the following cases for when to connect $`(z_1, b_1)`$ and
$`(z_2, b_2)`$:

a) if $`h(b_1,b_2) = 0`$ and $`g(z_1,z_2) = 1`$
b) if $`h(b_1,b_2) = 1`$
c) if $`h(b_1,b_2) = 2`$
d) if $`h(b_1,b_2) = 3`$ and $`g(z_1,z_2) \ne 1`$.

You can show this by again adding pairs of coordinate 2’s, one of
which is in the first position.

Note that this means that the $z$ coordinate, which is a member of
$`\mathbb{Z}/4`$ uses the Lee metric as distance.  One can use the
*Gray Map* 0->00, 1->01, 2->11, 3->10, which will convert the Lee
metric into the Hamming metric.  Thus the graph that we're interested
in has vertices labeled with elements of $`\mathbb{F}_2^n`$.  For
convenience define the map $`\phi: \mathbb{F}_2^{n+1} \rightarrow
\mathbb{F}_2^{n-1}`$, and $`\psi: \mathbb{F}_2^{n+1} \rightarrow
\mathbb{F}_2^2`$ as projection to the last $n-1$ coordinates, and
first 2 coordinates, respectively.  If we define the set $S$ as being
the union of the following sets

$`\{z: \psi(z) \in \{(0,1),(1,0)\}, \phi(z) = 0\}`$,
$`\{z : \text{wt}(\phi(z)) = 1,2\}`$
$`\{z: \psi(z) \in \{(0,0),(1,1)\}, \text{wt}(\phi(z)) = 3\}`$,

then edge in the graph exists between $z$ and $w$ if and only if
$z+w \in S$.
A large group of automorphisms of the graph is generated by
all permutations of the last $n-1$ coordinates, and all permutations
of the first 2 coordinates, along with adding an arbitrary member of
$`\mathbb{F}_2^{n+1}`$.  Thus the graph is vertex transitive.  Without
loss of generality, up to action by this group, we may assume that a
maximum independent set contains $0$, and does not contain $S$.  Note
that the vertex $(1,1,0, \dots, 0)$ has the same neighbors as 0.  Thus
we may remove both 0 and that vertex plus their common neighbors from
the graph.  For the resulting graph we have an automorphism group of
size $4 (n-1)!$, which consists of arbitrary permutations of the last
$n-1$ coordinates, permuations of the first 2 coordinates, plus adding
$(1,1)$ to the first two coordinates.  This group leaves the set
$\{(0, \dots, 0), (1,1,0, \dots, 0)\}$ invariant.

