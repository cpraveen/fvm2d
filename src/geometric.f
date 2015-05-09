C Calculate geometric quantities like face areas, normals, volumes, etc.
      subroutine geometric(coord, elarea, cvarea, ds, dsb, mcarea, 
     +                     wd, drmin)
      implicit none
      include 'param.h'
      include 'gdata.h'

      integer          esup1(mesup*npmax), esup2(npmax+1),
     &                 psup1(mpsup*npmax), psup2(npmax+1)
      double precision coord(2, npmax), elarea(ntmax), cvarea(npmax), 
     &                 ds(2,nemax), dsb(2,npmax), drmin(npmax), 
     &                 tcoord(2,ntmax), mcarea(npmax), wd(npmax)

c     Make sure elements are oriented ccw
      call tri_orient(elem, coord)

c     Find elements surrounding a point
      call el_surr_point(elem, esup1, esup2)

c     Find points surrounding a point
      call pt_surr_pt(esup1, esup2, elem, psup1, psup2)

c     Create edges
      call create_edge(psup1, psup2, edge)

c     Find element adjoining each edge
      call el_surr_edge(esup1, esup2, elem, edge, edneigh)

c     Write grid in gnuplot format for visualization
      call write_grid(coord, edge, edneigh)

c     Compute tcoord: centroid for Median cell, circuncenter for barth cell
      call tri_coord(coord, tcoord, elem)
 
c     Calculate element and control volume areas
      call areas(coord, tcoord, elem, elarea, cvarea, mcarea)

c     Find wall distance for Spalart-Allmaras model
      if(iflow .eq. turbulent) call wall_dist(edge, ptype, coord, spts,
     &                                        psup1, psup2, wd)

c     Find length of control volume faces
      call flength(coord, tcoord, edge, edneigh, ds, dsb)

      call reorder_edges(iedge, ds)

c     From this point, some edge operations may require iedge
c     Note that ne may not equal netot, be careful
c     If you move any of the subroutines be careful to use ne/netot

c     Length scale for time-step calculation
      call dtlength(coord, elarea, elem, drmin)

c     Identify boundary edges and make an index of them in beindx
      call make_bdedges(ptype, edge, edneigh, beindx)

c     Write dual grid into DUAL.DAT, use gnuplot to visualize it
      call write_dual(coord, tcoord, edge, iedge, edneigh)

      call geom_stat(elarea, cvarea, mcarea, ds)

      return
      end
