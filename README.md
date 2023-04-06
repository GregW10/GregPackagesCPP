# GregPackagesCPP
This is the repository which contains all of my frameworks (functions and classes) to simulate N-body systems in 3D space.

Every single line of code is 100% my own work. In addition, no third-party libraries/frameworks have been used anywhere and, therefore, everything is written in standard C++20. The code is cross-platform and will compile on Mac, Linux and Windows.

All files not beginning with "greg" are only helper files and are not part of the actual framework (they can be ignored). Nonetheless, they make use of my framework files and, thus, provide good examples of how my numerous classes and functions can be used.

## gregalg.hpp

This file is relatively small. I use it to place plenty of my miscellaneous algorithms that do not have a place in other header files.

Importantly, though, the file contains the definition of my `gtd::isNumWrapper` concept.

This concept encompasses any data type that behaves as if it were a number/fundamental type (i.e., it has `+`, `-`, `/`, `*` overloaded, among various others).

The idea behind `gtd::isNumWrapper` is that my templated N-Body classes (see below) should be able to use any data type that has the necessary overloads defined. The motivation for this is to allow an arbitrary-precision floating point type to be used within the N-Body classes, to allow for higher accuracy in the evolution of N-Body systems.

Nonetheless, support for this in my classes still needs work, given that the C/C++ standard library mathematical functions can only operate on fundamental types, which is a limiting factor in all calculations (e.g., `sqrtl` is used everywhere). Where these functions are used, I have always used the `long double` versions to maximise precision, but this is not a solution.

In the future, I will provide some mechanism of using custom functions that operate on the given data type to be used, which will then preserve the precision of any arbitrary precision type to be used within my N-Body classes.

Side note: I plan on writing my own arbitrary precision class, and will probably name the file `gregnum.hpp`.

## gregmat.hpp

This header file contains the `gtd::matrix<gtd::isNumWrapper T>` class.

`gtd::matrix`, as the name suggests, is a class whose objects represent mathematic matrices.

Appropriate overloads have been defined, such as addition `+` and multiplication `*`, as well as numerous useful methods, such as `gtd::matrix<T>::transpose`, `gtd::matrix<T>::determinant` (beware of stack overflows for this one as it uses recursion), `gtd::matrix<T>::has_inverse`, `gtd::matrix<T>::make_identity`, among others.

A number of static methods also exist, such as `gtd::matrix<T>::get_3D_rotation_matrix`, which returns a matrix used for rotating a 3D vector about one of the three coordinate axes.

## gregvec.hpp

This header file contains the `gtd::vector<gtd::isNumWrapper T>` abstract base class (I have marked this abstract class for removal, however, as it serves no purpose), the `gtd::vector2D<gtd::isNumWrapper T>` class and the `gtd::vector3D<gtd::isNumWrapper T>` class.

`gtd::vector2D<>` and `gtd::vector3D<>` objects represent 2D and 3D vectors in space, respectively.

I truly have gone to town with the number of methods and operator overloads I have defined for these classes, allowing a vast number of different vector operations to be performed, such as `gtd::vector3D<>::rodrigues_rotate` (for rotating the vector object about another vector by a given angle), `gtd::vector3D<>::magnitude`, `gtd::vector3D<>::unit_vector`, `gtd::vector3D<>::apply` (which applies a 3x3 `gtd::matrix<>` object to the vector), `gtd::vector<>::set_length`, among others.

All expected operator overloads are present, such as `+`, `-`, and `*` (the latter is for the dot product), as well as others which are unphysical (such as `/` between two vectors), but have been included for ease of use of the class.

There is a nested namespace with `gtd` called `gtd::vec_ops`, which, thus far, contains the `gtd::vec_ops::distance`, `gtd::vec_ops::cross` and `gtd::vec_ops::angle_between` functions.

## gregbod.hpp

This file contains the `gtd::body<isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t recFreq>` class used in all N-Body simulations.

A `gtd::body<>` object is one which represents a body in 3D space - it has a current position, velocity and kinetic energy (among others). It also has a given mass and radius.

The `M` template parameter represents the data type used to store the mass of a `gtd::body<>` object, `R` for the radius and `T` for the components of position and velocity - the position and velocity are stored as `gtd::vector3D<T>` objects.

A `gtd::body<>` object has the ability to store its "history" of previous positions, velocities and kinetic energies, with a frequency specified by the `recFreq` non-type template parameter (i.e., a `recFreq` of 3 means that every 3 updates of position and velocity will result in a recording).

`recFreq` was made a non-type template parameter to allow a template specialisation for `recFreq == 0` (i.e., no recording of history). This specialisation results in a much smaller static size of `gtd::body<>`, given that all the variables and data structures (`std::vector<>`s) used for storage can be removed in the specialisation.

### Note

This header file also contains the abstract `gtd::body_counter` (no template parameters) class, which is the parent class of `gtd::body<>`.

The purpose of the `gtd::body_counter` class is to ensure each `gtd::body<>` object is assigned a unique ID, such that no two `gtd::body<>` objects in the entire running of a program can be assigned the same unique ID. `gtd::body_counter` was purposefully made not to include any template types, so that it is a parent of **ALL** `gtd::body<>` classes (i.e., all versions of the `gtd::body<>` class - with different template parameters - inherit from the same `gtd::body_counter` class). This ensures that all instances of different `gtd::body<>` types still have unique IDs.

## gregsys.hpp

This is a biggie.

The `gregsys.hpp` header file contains the `gtd::system<gtd::isNumWrapper M, gtd::isNumWrapper R, gtd::isNumWrapper T, bool prog, bool mergeOverlappingBodies, int collisions, uint64_t memFreq, uint64_t fileFreq, bool binaryFile>` class.

A `gtd::system<>` object represents a N-Body system in 3D space, composed of numerous `gtd::body<>` objects.

The `M`, `R` and `T` template types represent the same as in `gtd::body<>`. `prog` represents whether the progress of evolution of the system should be printed to `stdout`. `mergeOverlappingBodies` represents whether overlapping bodies should cause an exception to be raised or should simply be merged either during the construction of a `gtd::system<>` object or when bodies are added - **this does *NOT* concern overlaps within the evolution itself**. `collisions` represents the collision resolution method to be employed within the `gtd::system<>` object. `memFreq` is the equivalent of `recFreq` in the `gtd::body<>` class (i.e., it represents the frequency of storage of position, velocity and kinetic energy values performed by the `gtd::body<>` objects themselves). `fileFreq` is similar to `memFreq`, but instead specifies the frequency of writing data to a file (not yet complete) - this is managed by the `gtd::system<>` object, however (and not by the `gtd::body<>` objects as in the case of `memFreq`). Finally, `binaryFile` represents whether the data being written to the file (if `fileFreq > 0`) is in binary format (as a `.nbod` file) or not.

The choice to make all of the above non-type template parameters (and not simply values passed to a constructor or methods) was to avoid any runtime overhead of repeatedly checking conditions and/or avoiding code bloat from producing loops that are almost identical but differ in a few lines, based on a parameter being true or false.

A `gtd::system<>` object works by passing it certain values upon construction (such as the time-step of evolution, the units to be used (which affects the numerical values calculated) and the number of iterations to perform). `gtd::body<>` objects are either passed to a `gtd::system<>` object in its constructor or later (via numerous methods). The `gtd::system<>::evolve` method carries out the evolution of all the bodies in the system via the integration method and force approximation technique passed to the method call (given as `static constexpr int` values that are part of the class).

After evolution, the positions and velocities of the bodies remain as they are, and if `recFreq != 0`, they can optionally be reverted to their state before the evolution.

If `recFreq != 0`, the `gtd::system<>::write_trajectories` method can be called to output the entire history of the system (body positions, velocities and energies, as well as total system kinetic energies, potential energies and total energies over time) to a .csv file (I will later include support for outputting to a binary `.nbod` file).

I have created a binary file format with the `.nsys` extension which represents a "snapshot" of a `gtd::system<>` object at a given time-step. This file format stores important properties of the `gtd::system<>` object itself (such as the time at which the snapshot was taken, the units for mass, distance and time, the time-step, etc...) as well as the positions and velocities of every single body in the system.

I have written a method which allows a `gtd::system<>` object to output a `.nsys` file of itself called `gtd::system<>::to_nsys`. To read in a `.nsys` file, a `gtd::system<>` object must be constructed from a `.nsys` file by passing its path to the constructor (certain non-framework helper files such as `scat.cpp` make use of this functionality).

The force on body $1$ due to body $2$ is calculated as:

$$\mathbf{F_{12}} = \frac{G m_1 m_2}{r_{12}^2} \mathbf{\hat{r}_{12}}$$

where $\mathbf{\hat{r}_{12}}$ points from body $1$ to body $2$.

Dividing by the mass of body $1$ results in the acceleration of body $1$ due to body $2$:

$$\mathbf{a_{12}} = \frac{Gm_2}{r_{12}^2}\mathbf{\hat{r}_{12}}$$

This can now be trivially generalised to calculate the total force on body $i$ due to the other $N - 1$ bodies in the system:

$$\mathbf{a_{i}} = G\sum_{j = 0,j\neq i}^{N - 1}{\frac{m_j}{r_{ij}^2}\mathbf{\hat{r}_{ij}}}$$

This results in a total of $\frac{N(N - 1)}{2}$ force calculations in direct integration methods (hence the $\mathcal{O}(N^2)$ time complexity).

Calculating the acceleration $\mathbf{a_i}$ of each body is used in time-stepping the system of bodies forward in time. For example, for the Euler integration method, the acceleration is calculated once for a given body at the start of a step, and then used to calculate its velocity at the next time-step, whilst its velocity is used to calculate its position at the next time-step:

$$\mathbf{v_{n+1}} = \mathbf{v_n} + \Delta t\mathbf{a_n}$$

$$\mathbf{r_{n+1}} = \mathbf{r_n} + \Delta t\mathbf{v_n}$$

where $\mathbf{r_n}$, $\mathbf{v_n}$ and $\mathbf{a_n}$ are the position, velocity and acceleration, respectively, of a body at step $n$, $\mathbf{r_{n+1}}$ and $\mathbf{v_{n+1}}$ are the position and velocity, respectively, of a body at step $n + 1$ (one time-step in the future), and $\Delta t$ is the time-step.

For the Kick-Drift-Kick (KDK) and Drift-Kick-Drift (DKD) Leapfrog integration methods, only a single force evaluation per time-step is also required, but these are 2<sup>nd</sup> order accurate (as opposed to 1<sup>st</sup> order accurate, as in the case of the Euler method). Other methods, such as the midpoint method and modified Euler method, require 2 force evaluations per time-step, and are also 2<sup>nd</sup> order accurate.

The above 5 integration methods are currently available in the `gtd::system<>` class. In the future, higher-order methods such as the 4<sup>th</sup> order Runge-Kutta method and the 4<sup>th</sup> order Hermite scheme will be employed. These offer higher accuracy, but also require increased computation time, given that they require more force evaluations per time-step.

Amongst other additions, to avoid singularities caused by two bodies overlapping perfectly (such that $r_{ij} = 0$) or extreme forces when $r_{ij}$ is very small, I will add a softening parameter $\epsilon$, which will alter the equation for the force between two bodies:

$$\mathbf{F_{12}} = \frac{G m_1 m_2}{r_{12}^2 + \epsilon^2} \mathbf{\hat{r}_{12}}$$

$\epsilon$ should be a small quantity, such that its effect is minimal for characteristic distances between bodies, but will dominate when $r_{ij} \ll \epsilon$ and prevent the force from tending towards infinity (which is unphysical).

## gregbmp.hpp

This file contains the parent class of `gtd::astro_scene<>` (which is the main rendering class in **gregastro.hpp**), namely `gtd::bmp`.

The `gtd::bmp` class provides an intuitive way of generating and manipulating BMP images.

Creating a .bmp image can be done in two lines with my class:

```
gtd::bmp bmp_image;
bmp_image.write();
```

Of course, the above would just produce a black `2400x1600` image, but the class contains numerous methods, such as `gtd::bmp::set_pixel` and `gtd::bmp::fill_rect` which allow manipulation of the .bmp image.

The image stores the internal pixel array on the heap (which it allocates upon construction), which it then writes to a file when `gtd::bmp::write` is called (the method generates a default file path using the current date & time if no path is passed).

There are also other classes, such as `gtd::pixel` and `gtd::square` which are used in the `gtd::bmp::draw_circle` method (which was a lot of work) to produce an anti-aliased circle. My circle anti-aliasing algorithm is slow and, thus, subject to future improvements, but works well. Basically, a rough, aliased circle is initially drawn, and then the area of the pixels on the edge of the circle that is covered by the circle is approximated for each pixel by sampling 256 points placed within each pixel. This area is then used to shade the pixel appropriately (the shading is directly proportional to the area covered). The shading also takes into account the background colour (surrounding the circle), so that the pixel on the edge of the circle ends up with a colour between the colour of the circle and the background.

### Note

The `gtd::bmp` class contains the raw pixel array (of `gtd::color`s) as a protected member named `data`. Thus, subclasses of `gtd::bmp` can directly access and manipulate the underlying pixel array, without having to worry about the I/O of the .bmp file at all.

`data` is a double pointer of type `gtd::color**`, where `gtd::color` is a record type (`struct`), implemented like so:

```
namespace gtd {
#pragma pack(push, 1)
  struct color {
    unsigned char b;
    unsigned char g;
    unsigned char r;
  };
#pragma pack(pop)
}
```

The `#pragma pack` preprocessor directive is vital in this context, as it removes the padding that would be added to `gtd::color` by default. Hence, `sizeof(gtd::color) == 3`. The colours are in BGR order and not RGB order, since .bmp files use the BGR format.

I have also created a namespace within `gtd` called `colors`. This namespace contains common colours, such `gtd::colors::aqua`, `gtd::colors::azure`, `gtd::colors::yellow`, `gtd::colors::pink`, `gtd::colors::forest_green`, and many more.

## gregastro.hpp

This file contains numerous classes that are responsible for rendering .bmp images of the simulations using ray-tracing. The main rendering class is `gtd::astro_scene<gtd::isNumWrapper M, gtd::isNumWrapper R, gtd::isNumWrapper T, gtd::isNumWrapper PosT, gtd::isNumWrapper DirT, gtd::isNumWrapper DistT, gtd::isNumWrapper LenT, gtd::isNumWrapper LumT, bool modulatedBrightness, uint64_t recFreq>` and is a child of `gtd::bmp`.

The basic premise is that a camera is treated as a "box" in 3D space, with a pinhole on one end and a receptor covering the opposite end that is made up of all the pixels composing the image. Line-of-sight (LOS) rays are "shot out" from each pixel on the receptor and **through the pinhole** towards the scene (this is taken care of by my `gtd::camera<gtd::isNumWrapper PosT, gtd::isNumWrapper DirT, gtd::isNumWrapper DisT>` class), and intersection with bodies in the scene is calculated (using the equation of a ray in 3D space and the equation for a sphere - this breaks down to solving the quadratic formula: 2 solutions indicates the ray enters through one side of the sphere and exits out the back; 1 solution indicates a ray is tangent to a sphere; 0 solutions indicates the ray does not intersect).

The equation for a ray is given as such:

$$\mathbf{r} = \mathbf{p} + l\mathbf{\hat{u}}$$

where $\mathbf{p}$ is its origin, $\mathbf{\hat{u}}$ is the direction it is pointing in and $l$ is its length, which combine to give $\mathbf{r}$, which is its end point. $l$ is what is solved for, given a number of spheres within a scene. The equation of a sphere is given by:

$$(x - \alpha)^2 + (y - \beta)^2 + (z - \gamma)^2 = R^2$$

where $\alpha$, $\beta$ and $\gamma$ are the $x$, $y$ and $z$-coordinates of the sphere's centre, respectively, and $R$ is its radius. Thus, in the equation, $x$, $y$ and $z$ denote the coordinates of any point on the sphere's surface. Therefore, given a sphere at coordinates $(\alpha, \beta\, \gamma)$, the possible intersection point of a ray with this sphere can be calculated by plugging the $x$, $y$ and $z$ components of the ray into $x$, $y$ and $z$ of the equation for the sphere, to end up with an equation containing only one unknown, $l$ (the length of the ray):

$$r_x = p_x + l\hat{u}_x$$

$$r_y = p_y + l\hat{u}_y$$

$$r_z = p_z + l\hat{u}_z$$

$$(p_x + l\hat{u}_x - \alpha)^2 + (p_y + l\hat{u}_y - \beta)^2 + (p_z + l\hat{u}_z - \gamma)^2 - R^2 = 0$$

The hat is kept above the components of the $\mathbf{\hat{u}}$ direction unit vector to emphasize that these are components of a unit vector.

This equation can then be rewritten as a quadratic equation:

$$al^2 + bl + c = 0$$

with $a$, $b$ and $c$:

$$a = \hat{u}^2 = 1$$

$$b = 2(\mathbf{p} \cdot \mathbf{\hat{u}} - \alpha\hat{u}_x - \beta\hat{u}_y - \gamma\hat{u}_z) = 2(\mathbf{p} - \mathbf{r_s})\cdot \mathbf{\hat{u}}$$

$$c = p^2 - 2(\alpha p_x + \beta p_y + \gamma p_z) + \alpha^2 + \beta^2 + \gamma^2 - R^2 = p^2 - 2\mathbf{r_s}\cdot \mathbf{p} + {r_s}^2 - R^2$$

where $\mathbf{r_s} = \alpha\mathbf{\hat{i}} + \beta\mathbf{\hat{j}} + \gamma\mathbf{\hat{k}}$ (position vector of the centre of the sphere), $\hat{u}^2 = \mathbf{\hat{u}}\cdot \mathbf{\hat{u}}$, $p^2 = \mathbf{p}\cdot \mathbf{p}$ and ${r_s}^2 = \mathbf{r_s}\cdot \mathbf{r_s}$ (squares of the magnitudes of the vectors).

Finally, these quantities can be plugged into the quadratic formula:

$$l = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

to solve for $l$. As mentioned above, 2 solutions indicates two intersections with the sphere, 1 indicates a tangent ray, and 0 means the ray does not intersect with the sphere (and continues on towards $\infty$).

The intersection of a ray shot out from a given pixel is calculated for all bodies in the scene: if the ray does not inersect with any (i.e., if there are 0 solutions to the quadratic equation for each body), then the pixel is coloured black, else, the body that is closest to the camera is selected as the intersected body (i.e., the one which results in the ray having the shortest length $l$) and the pixel is coloured as the colour of the intersected body (either as set by the user or given by the randomly generated colour that was assigned to the body when it was added to the `gtd::astro_scene<>` object).

Next, the brightness of the pixel must be calculated (since the above just determined the intersected body and, hence, the colour). On top of the raw array of `gtd::color`s (accessed via the `data` double pointer of the parent `gtd::bmp` class), there is a **second array of `LumT` values** with the same dimensions as the colors array, which **stores the brightness of each pixel**.

Light is "produced" inside a scene via `gtd::star<gtd::isNumWrapper M, gtd::isNumWrapper R, gtd::isNumWrapper T, gtd::isNumWrapper PosT, gtd::isNumWrapper DirT, gtd::isNumWrapper LenT, gtd::isNumWrapper LumT, uint64_t recFreq>` objects. `gtd::star<M, R, T, PosT, DirT, LenT, LumT, recFreq>` is a subclass of `gtd::body<M, R, T, recFreq>`. A `gtd::star<>` object contains 1 or more `gtd::light_src<gtd::isNumWrapper PosT, gtd::isNumWrapper DirT, gtd::isNumWrapper LenT, gtd::isNumWrapper LumT>` objects (if it contains 1 `gtd::light_src<>`, this source is treated as being at the centre of the star, else, the sources are distributed randomly across the surface of the `gtd::star<>`, producing more realistic shadows for large numbers of light sources). A `gtd::star<>` object contains a given luminosity (whose type is given by the `LumT` template parameter) that is distributed equally among all its `gtd::light_src<>` objects (i.e., for $N$ `gtd::light_src<>` objects, each source has a luminosity of $\frac{L_{star}}{N}$).

The brightness of each pixel (for a pixel that has produced a LOS ray which intersected with a body) is determined via these steps:

1. For each `gtd::light_src<>` object in the scene, cast a light ray from the light source to the point at which the LOS ray intersected with a body.
2. If the `gtd::light_src<>` object resides on a `gtd::star<>` with more than 1 light source, check that the star does not block the ray by taking the dot product of the ray's direction with the normal to the star at the location of the light source - it must be **positive**, else, continue on to the next `gtd::light_src<>`.
3. Check that the (closest) `gtd::body<>` the light ray intersects with is the same as the one the LOS ray intersects with, else, move on to the next `gtd::light_src<>`.
4. Check that the point the LOS ray intersects with is the "entrance" point of the light ray on the sphere, and not the "exit" point (in which case the point would not be illuminated by that light source) by taking the dot product of the normal to the sphere at that point with the light ray and ensuring it is **negative** - if not, move on to the next `gtd::light_src<>`.
5. By this stage it is certain that the light ray from the given `gtd::light_src<>` intersects with the point on the sphere intersected by the LOS ray emitted by the `gtd::camera<>` object, so the brightness (flux) at the point due to the given light source is calculated using the inverse square law and then multiplied by the negative dot product of the normal at the point on the sphere with the direction of the light ray (since the dot product itself is negative) - this was already calculated in step 4.
6. Finally, add this "flux" or "brightness" onto the cumulative flux value stored in the 2D array of brightness values (for the corresponding pixel).

The normal to sphere at the given point is calculated by:

$$\mathbf{n} = \mathbf{p_{lr}} + l_{lr}\mathbf{\hat{u}_{lr}} - \mathbf{r_s}$$

where $\mathbf{p_{lr}}$, $l_{lr}$, and $\mathbf{\hat{u}_{lr}}$ are the origin, length, and direction of the given light ray, respectively, and $\mathbf{r_s}$ is the position vector of the centre of the sphere.

Thus, the brightness at the point due to a single light source is calculated as:

$$I = \left( \frac{L}{l_{lr}^2} \right) \left| \mathbf{\hat{n}}\cdot \mathbf{\hat{u}_{lr}} \right|$$

where $\mathbf{\hat{n}}$ is the unit vector of $\mathbf{n}$.

Although flux (without the adjustment due to the angle of incidence) is technically given by:

$$f = \frac{L}{4\pi d^2}$$

(where $d$ is the distance to a light source), since all brightnesses are normalised in `gtd::astro_scene<>` (to ensure the maximum R, G, or B value is 255 - thus not reducing the limited 0-255 range any further) the $4\pi$ is dropped and the proportionality is used instead to improve efficiency:

$$f \propto \frac{L}{d^2}$$

Once the steps above have been taken for all pixels, all brightness values in the array are normalised by dividing by the maximum brightness value - such that all brightness values end up between $0$ and $1$. Finally these normalised brightness values are multiplied by the B, G, and R channels of the `gtd::color`s in the main colour array to produce the final image with the correct brightness modulations.

However, there is a small subtlety to finding the maximum brightness value. I have made the entire ray-tracing and rendering algorithm multithreaded, with the number of threads determined using `std::thread::hardware_concurrency`. I split the image up initially by rows (horizontally) - since color and brightness values for pixels are stored contiguously by rows and not by columns (thus minimising the amount of "jumps" in memory that have to be performed) - and only start to split the image up by columns (vertically) if the number of threads is greater than the number of rows in the image (highly unlikely for an average computer and average-sized image). Thus, all threads end up with an area of the image to tackle that is almost equal (and precisely equal if, in the first case, the number of threads is a factor of the height of the image and, in the second case, if the height minus the number of threads is a factor of the height and if the number of threads divided by the height is a factor of the width).

Thus, I make all threads calculate the maximum brightness in their area of the image, then make them wait for each other to finish (using `std::latch`) and then add their maximum brightnesses to a `std::vector<LumT>` object that is accessed using a `std::mutex` (to avoid a race condition). Then, each thread notifies a completely separate thread (that is only concerned with calculating the maximum brightness from all the threads) - using a `std::condition_variable` object - to calculate the maximum brightness, and this separate thread only proceeds once the number of elements in the `std::vector<LumT>` object containing the threads' maximum brightnesses equals the number of threads. Once this separate thread has calculated the maximum brightness out of all the threads, the other threads are notified and use this maximum to normalise each brightness value in the array of brightnesses as they multiply them by the colours in the `gtd::color`s array.

I will attempt to add support for GPUs in the future, as `std::thread` is CPU-specific.

## greg8tree.hpp

This header file contains my Barnes-Hut classes.

The Barnes-Hut algorithm is one which manages to reduce the time complexity of a direct N-Body integration algorithm from $\mathcal{O}(N^2)$ to $\mathcal{O}(Nlog(N))$ by approximating the force from distant bodies on a given body as if the distant bodies were a single, larger body at their centre-of-mass (COM).

This is done by first selecting a bounding cube (which is large enough to fit all bodies), and then feeding bodies into this cube, one-by-one. The moment 2 bodies are found within a single cube, the cube is recursively divided into 8 sub-cubes until each body is within its own cube (so this first division happens the moment the second body is added to the root cube).

Thus, one ends up with an octree data structure in which each leaf node contains a single body (leaf nodes, in my version of the algorithm, are only created if a body will reside inside it, as such, there are no empty leaf nodes).

A salient feature of my `gtd::bh_cube<gtd::isNumWrapper M, gtd::isNumWrapper R, gtd::isNumWrapper T, uint64_t rF>` class (the main class representing a Barnes-Hut cube or node) is that it does **not** use recursive functions (apart from the `<<` overload for printing to `stdout`). Instead, iteration is used to perform all recursive operations. This greatly increases performance (as the overhead of setting up a new stack frame and creating local copies of variables at each new recursion depth is avoided) and reduces the chance of a stack overflow (although the latter was unlikely anyway, given the fact that one attraction of the Barnes-Hut algorithm is that the tree height never gets too large).

The `gtd::bh_cube<>` class contains three non-abstract nested iterator classes, namely `gtd::bh_cube<>::iterator` (for iterating over all bodies contained within a given `gtd::bh_cube` object), `gtd::bh_cube<>::nn_iterator` (for iterating over the nearest neighbours of a given body in a given cell) and `gtd::bh_cube<>::cell_iterator` (for iterating over all cells that satisfy the opening angle parameter, given a target body for which force will be calculated).

Thus far, my Barnes-Hut classes have not been fully implemented within my `gtd::system<>` class, but this will be amended in the near future.

## gregstr.hpp

This file is not directly relevant to N-Body simulations, but has proven to be highly useful, nonetheless.

In it, I have defined numerous useful functions, such as `gtd::get_date_and_time`, `gtd::get_home_path` and `gtd::strsplit`, among many others.

Most importantly, however, the file contains my `gtd::String` class. This class is similar to `std::string` in that it provides a user-friendly interface for dealing with strings, but has important differences.

As a side note, it is important to mention that, although similar to `std::string`, `gtd::String` does not make use of a `std::string` anywhere to provide functionality, nor is it implemented in the way `std::string` is (I presume, since I have never actually looked at the `std::string` implementation, although I know how it works).

The main difference is that `gtd::String` stores the characters representing the string in the **middle** of its `char` array (and not at the beginning, as with `std::string`). This allows $\mathcal{O}(1)$ (constant) time complexity for addition of characters at the end **and** at the beginning of the string (unless a reallocation has to be performed).

Furthermore, `gtd::String` provides much functionality that `std::string` does not, whilst still providing its own implementation of almost all functions in `std::string`, thus allowing a familiar user-experience.