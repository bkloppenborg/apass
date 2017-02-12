Hi Arne,

Good news. I didn't think about the merging fred file problem at all this weekend until around 8 PM tonight at which point I sat down and clunked out the solution in 11 pages of notes and pseudcode. The solution is highly parallel, extremely fast, and will keep residual data around for subsequent changes. The main steps are A, AA (to keep consistency with my notes), B, C, and D.

Please give it a quick read to see if I've missed something obvious. I think the only issue I might have is dealing with the poles, but I might project that data into a Cartesian coordinate system with the pole as the origin.

Regards,

Brian

-----

A. Subdivide the sphere into zones (as you suggested) and insert data into zone file

I'll need your advice here, but I think we could shoot for something like 5x2.5 degree zones except at the poles which would have one merged zone for each pole. We'll have to figure out the size for that as we need to keep the zones small enough to fit into memory. The general method is as follows:

    Create a quadtree which divides the whole sphere to depth L. Each leaf node will get a unique ID which will correspond to the zonefile (a .fredbin). Reserve id=0,1 for the polar nodes.
    Modify the quadtree such that nodes with dec < X or dec > X points to the polar nodes. Save a copy of the quadtree to disk.
    Initialize a .fredbin (binary fred) file for each zone.

AA: Insert FRED data

    For each fred file, insert the data into a zone file. Due to the spatial locality of the data, open zonefiles as needed and close them once a fredfile is done or the maximum number of open zonefiles is reached.

Input: fred files

Output: A global quadtree storing zone data. Several zone files (.fredbin)

Other: After steps 1-3 are done, the operations can be done in parallel provided that we write some lockfiles to protect the .fredbin files. This process will be limited by disk I/O and memory.

-----

B. Generate rectangles in each zonefile

I've written some simple code that does this step using brute force searching (e.g. O(N^2)). The quadtree should give me O(N log_4 N) characteristics.

    1. Create a quadtree for the zone. The quadtree will subdivide the zone to some depth M. I'm thinking 10", but something bigger might be ok. I call the leaves of this quadtree a "cell".
    For each star in the zonefile (this is an algorithm I have written down in pseudocode)
        create a bounding rectangle with 1" on each side
        traverse the quadtree to a cell
        gather a list of overlapping rectangles
        remove the overlapping rectangles from the cell
        merge the list of overlapping rectangles
        append the resulting rectangle to the cell
    Merge rectangles within the zonefile
        foreach rectangle in a quadtree cell
            if the rectangle is on the edge of the zone
                save its ID to a "zoneedge" file
            if it is within radius R of the edge (we'll need to define this)
                find the adjacent cells (I have this algorithm written out... it's kinda lengthy, but easy to implement)
                foreach cell
                    if any of the rectangles overlap,
                        merge data into this cell's rect. Or we could follow some other convention (like lesser RA, more positive DEC)
                        mark other rects as moved to (zone_id, rect_id). This leaves a "pointer" in the tree for future operations.
            save the rectangle's data to a zonefile-rectid.fredbin file

Input: Zone files (.fredbin)

Output: A zone quadtree storing a list of spatially organized rectangles. Several zone_id-rectangle_id.fredbin files. A list of rectangles on the edge of the zone (a "zoneedge" file).

Other: This can all be done in parallel as everything is zone-independent! This process will be limited by disk I/O and system memory.

-----

C. Merge Rectangles on the edge of zones

    load the global quadtree
    foreach rect in a zoneedge file
        find adjacent zones (uses the same three traversal method as above except we handle the RA = 0 -> 360 and polar wraparound
        load zone quadtree
        find adjacent cells
        foreach rect in cell
            if there is an overlap
                merge data into this cell's rect (or follow the other convention)
                mark moved to (zone_id, rect_id). Again, this leaves a "pointer" in the tree for future operations.
                write the changes to disk
        save and close the zone quadtrees (save memory, these things are probably big)

Input: zoneedge files

Output: zone_id-rectangle_id.fredbin files

Other: This is a serial process.

-----

D. Statistics!

At this point we have a bunch of these zone_id-rect_id.fredbin files on disk. If all worked ok, it will contain the data for one (and only one) star. There might be some exceptions, for which we'll have to check, but that will be really easy. Now we simply process these files using standard statistical methods and whatever filtering based upon CCD location, sky conditions, etc.

-----

E. Other

I've also written down algorithims for adding more data and removing data. These just involve loading the global quadtree, traversing to a zone, loading the zone's quadtree, traversing to a cell, and then looking up the corresponding rectangle's data.

We'll have to figure out how to mangle the results back into a format your subsequent code can understand, but that too should be easy.


