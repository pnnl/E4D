module buildmesh2

  use vars
  use input

  implicit none

  private
  public :: build_tetgen_direct

contains

!__________________________________________________________________________________________________
  subroutine build_tetgen_direct(i_flag)

    !! Build a TetGen PLC directly from the mesh configuration file.
    !!
    !! This routine removes the Triangle -> surface.1.node -> TetGen handoff that the
    !! original build path used.  Instead, it:
    !!   1) reads the mesh cfg file,
    !!   2) reorders the control points so that surface points come first,
    !!   3) triangulates the top surface directly in memory,
    !!   4) writes one TetGen .poly file, and
    !!   5) calls TetGen on that file.
    !!
    !! The goal is to make the surface/volume interface explicit in one place so the
    !! mesh build is easier to reason about and less dependent on post-processing.

    implicit none

    logical, intent(in) :: i_flag

    integer :: i, j, k, m, ios, npre, nchar
    integer :: n_points, nsurf, nedge, nint, nplc, nholes, nzones
    integer :: tetflag, exoflag, exoflag2d, writit
    integer :: nsurf_seeds, boundary_seed_count
    integer :: nbnd_facets, nside_facets, ntop_facets, nplc_facets, nfacets
    integer :: d1, d2, d3, el1, el2, el3, el4, el5
    integer :: nedge0, ecount, tcount, nsp2
    integer :: icount, idx, top_id, bottom_id, facet_count
    integer :: a, b, c, p
    integer :: ntop, ntri
    integer :: ntri_max
    integer :: tri_i, tri_j, tri_k
    integer :: n_fallback
    integer :: nnodes
    integer :: local_id
    integer :: surf_hit
    integer :: ncmax

    integer, dimension(:), allocatable :: tflag, edge_flag, pmap
    integer, dimension(:), allocatable :: edge_seq, boundary_order
    integer, dimension(:), allocatable :: surf_index, interior_index
    integer, dimension(:), allocatable :: tri_v1, tri_v2, tri_v3
    integer, dimension(:), allocatable :: tri_work
    integer, dimension(:), allocatable :: plc_npts, plc_bnd
    integer, dimension(:,:), allocatable :: plc_pts
    integer, dimension(:,:), allocatable :: zone_labs

    real*8 :: xorig, yorig, zorig, min_elev, miny_ap, maxy_ap
    real*8 :: x1, y1, x2, y2, x3, y3, area2, tol
    real*8 :: dist, w1, w2, w3
    real*8, dimension(:,:), allocatable :: tall_points, all_points
    real*8, dimension(:,:), allocatable :: holes_xyz, zones_xyz
    real*8, dimension(:), allocatable :: edge_sep, hole_vols, zone_sig, zone_sigi
    real*8, dimension(:), allocatable :: surface_x, surface_y, surface_z
    real*8, dimension(:), allocatable :: tri_x1, tri_y1, tri_x2, tri_y2, tri_x3, tri_y3
    real*8, dimension(:), allocatable :: point_x, point_y

    character*20 :: qual, max_vol
    character*80 :: tetcom, exocom, tricom, exocom2d
    character*40 :: mesh_cfg, mesh_poly, mesh_prefix
    logical :: exst, mtran
    logical :: found_ear, inserted, boundary_point, surface_point

    open(52, file='mesh_build.log', action='write', status='replace')
    close(52)

    !! The mesh file is supplied through the global variables module.
    inquire(file=trim(mshfile), exist=exst)
    if (.not. exst) stop 'Mesh configuration file not found.'

    mesh_cfg = trim(mshfile)
    mesh_poly = trim(mshfile)
    mesh_prefix = trim(mshfile)

    !! Read the mesh configuration file.
    open(10, file=mesh_cfg, status='old', action='read')
    read(10, *, iostat=ios) qual, max_vol
    if (ios /= 0) stop 'Error reading quality/max volume line in mesh cfg.'
    read(10, *, iostat=ios) min_elev
    if (ios /= 0) stop 'Error reading minimum elevation in mesh cfg.'
    read(10, *, iostat=ios) tetflag
    if (ios /= 0) stop 'Error reading tetgen build flag in mesh cfg.'
    read(10, *, iostat=ios) tetcom
    if (ios /= 0) stop 'Error reading tetgen executable name in mesh cfg.'
    read(10, *, iostat=ios) tricom
    if (ios /= 0) stop 'Error reading triangle executable name in mesh cfg.'
    read(10, *, iostat=ios) n_points
    if (ios /= 0) stop 'Error reading number of control points in mesh cfg.'

    if (n_points < 3) stop 'Mesh cfg must contain at least three control points.'

    allocate(tall_points(n_points, 3), tflag(n_points), edge_flag(n_points), pmap(n_points))
    allocate(all_points(n_points, 3))

    !! Read the raw control points.  The cfg format uses a per-point flag:
    !!   1 = additional surface control point
    !!   2 = outer boundary point on the top surface
    !!   everything else = interior point used for refinement.
    tflag = 0
    edge_flag = 0
    pmap = 0
    do i = 1, n_points
       read(10, *, iostat=ios) j, tall_points(i,1:3), tflag(i)
       if (ios /= 0) stop 'Error reading control point block in mesh cfg.'
       if (tflag(i) == 2) edge_flag(i) = 1
    end do

    !! Reorder the points so that surface points come first.
    !! This keeps the top surface and the refinement points separate.
    nsurf = 0
    nedge = 0
    nint = 0
    do i = 1, n_points
       if (tflag(i) == 1) then
          nsurf = nsurf + 1
          all_points(nsurf,1:3) = tall_points(i,1:3)
          pmap(i) = nsurf
       end if
    end do
    do i = 1, n_points
       if (tflag(i) == 2) then
          nsurf = nsurf + 1
          nedge = nedge + 1
          all_points(nsurf,1:3) = tall_points(i,1:3)
          pmap(i) = nsurf
       end if
    end do
    do i = 1, n_points
       if (tflag(i) /= 1 .and. tflag(i) /= 2) then
          nint = nint + 1
          all_points(nsurf + nint,1:3) = tall_points(i,1:3)
          pmap(i) = nsurf + nint
       end if
    end do

    !! A simple mesh sanity check: a top surface with fewer than three boundary points
    !! cannot be triangulated.
    if (nedge < 3) stop 'At least three boundary points are required.'

    if (nsurf > 0) then
       allocate(surf_index(nsurf))
       do i = 1, nsurf
          surf_index(i) = i
       end do
    end if
    allocate(boundary_order(nedge))
    do i = 1, nedge
       boundary_order(i) = 0
    end do
    if (nint > 0) then
       allocate(interior_index(nint))
       do i = 1, nint
          interior_index(i) = nsurf + i
       end do
    end if

    !! Build the list of boundary points in order.
    !! The original code used a nearest-neighbor walk; we keep that idea here because
    !! the cfg files already provide the boundary as a simple closed loop.
    allocate(edge_seq(nedge), edge_sep(nedge))
    ecount = 0
    do i = 1, nsurf
       if (tflag(i) == 2) then
          ecount = ecount + 1
          boundary_order(ecount) = i
       end if
    end do

    edge_seq(1) = 1
    do i = 1, nedge - 1
       x1 = all_points(boundary_order(edge_seq(i)),1)
       y1 = all_points(boundary_order(edge_seq(i)),2)
       do j = 1, nedge
          edge_sep(j) = sqrt((all_points(boundary_order(j),1) - x1)**2 + &
                             (all_points(boundary_order(j),2) - y1)**2)
       end do
       edge_sep(edge_seq(1:i)) = maxval(edge_sep)
       edge_seq(i+1) = minloc(edge_sep(:),1)
    end do

    !! Convert the boundary order into surface-point indices.
    do i = 1, nedge
       edge_seq(i) = boundary_order(edge_seq(i))
    end do

    !! Read optional PLC, hole, region, and export options.
    read(10, *, iostat=ios) nplc
    if (ios /= 0) stop 'Error reading PLC count in mesh cfg.'

    if (nplc > 0) then
       allocate(plc_npts(nplc), plc_bnd(nplc))
       allocate(plc_pts(nplc, 256))
       plc_pts = 0
       do i = 1, nplc
          read(10, *, iostat=ios) plc_npts(i), plc_bnd(i)
          if (ios /= 0) stop 'Error reading PLC header in mesh cfg.'
          if (plc_npts(i) < 3) stop 'Each PLC needs at least three points.'
          if (plc_npts(i) > size(plc_pts,2)) stop 'PLC point list exceeds the fixed scratch buffer.'
          read(10, *, iostat=ios) plc_pts(i,1:plc_npts(i))
          if (ios /= 0) stop 'Error reading PLC point list in mesh cfg.'
          do j = 1, plc_npts(i)
             if (plc_pts(i,j) < 1 .or. plc_pts(i,j) > n_points) stop 'PLC point index out of range.'
          end do
       end do
    end if

    read(10, *, iostat=ios) nholes
    if (ios /= 0) stop 'Error reading hole count in mesh cfg.'
    if (nholes > 0) then
       allocate(holes_xyz(nholes,3))
       do i = 1, nholes
          read(10, *, iostat=ios) j, holes_xyz(i,1:3)
          if (ios /= 0) stop 'Error reading hole location in mesh cfg.'
       end do
    end if

    read(10, *, iostat=ios) nzones
    if (ios /= 0) stop 'Error reading zone count in mesh cfg.'
    if (nzones > 0) then
       allocate(zones_xyz(nzones,3), zone_labs(nzones,2), hole_vols(nzones), zone_sig(nzones))
       if (i_flag) allocate(zone_sigi(nzones))
       do i = 1, nzones
          if (i_flag) then
             read(10, *, iostat=ios) zone_labs(i,1), zones_xyz(i,1:3), hole_vols(i), zone_sig(i), zone_sigi(i)
          else
             read(10, *, iostat=ios) zone_labs(i,1), zones_xyz(i,1:3), hole_vols(i), zone_sig(i)
          end if
          if (ios /= 0) stop 'Error reading zone definition in mesh cfg.'
          zone_labs(i,2) = zone_labs(i,1)
       end do
    end if

    read(10, *, iostat=ios) exoflag
    if (ios /= 0) stop 'Error reading Exodus flag in mesh cfg.'
    read(10, *, iostat=ios) exocom
    if (ios /= 0) stop 'Error reading Exodus command in mesh cfg.'
    read(10, *, iostat=ios) i
    if (ios /= 0) stop 'Error reading mesh-translation flag in mesh cfg.'
    if (i == 0) then
       mtran = .false.
    elseif (i == 1) then
       mtran = .true.
    else
       stop 'Mesh translation flag must be 0 or 1.'
    end if
    close(10)

    !! If requested, translate the geometry to improve floating-point conditioning.
    !! The translation is applied after the point reordering so that all derived data
    !! share the same coordinate frame.
    xorig = 0.0d0
    yorig = 0.0d0
    zorig = 0.0d0
    if (mtran) then
       nsurf_seeds = 0
       do i = 1, nsurf
          if (tflag(i) /= 2) then
             xorig = xorig + all_points(i,1)
             yorig = yorig + all_points(i,2)
             zorig = zorig + all_points(i,3)
             nsurf_seeds = nsurf_seeds + 1
          end if
       end do
       if (nsurf_seeds > 0) then
          xorig = xorig / dble(nsurf_seeds)
          yorig = yorig / dble(nsurf_seeds)
          zorig = zorig / dble(nsurf_seeds)
       else
          xorig = sum(all_points(1:nsurf,1)) / dble(nsurf)
          yorig = sum(all_points(1:nsurf,2)) / dble(nsurf)
          zorig = sum(all_points(1:nsurf,3)) / dble(nsurf)
       end if
    end if

    do i = 1, n_points
       all_points(i,1) = all_points(i,1) - xorig
       all_points(i,2) = all_points(i,2) - yorig
       all_points(i,3) = all_points(i,3) - zorig
    end do
    min_elev = min_elev - zorig
    if (nholes > 0) then
       holes_xyz(:,1) = holes_xyz(:,1) - xorig
       holes_xyz(:,2) = holes_xyz(:,2) - yorig
       holes_xyz(:,3) = holes_xyz(:,3) - zorig
    end if
    if (nzones > 0) then
       zones_xyz(:,1) = zones_xyz(:,1) - xorig
       zones_xyz(:,2) = zones_xyz(:,2) - yorig
       zones_xyz(:,3) = zones_xyz(:,3) - zorig
    end if

    !! Save the translation vector so later steps can reconstruct the original frame.
    nchar = 0
    do i = 1, 40
       if (mesh_cfg(i:i) == '.') then
          nchar = i
          exit
       end if
    end do
    if (nchar == 0) nchar = len_trim(mesh_cfg)
    open(22, file=mesh_cfg(1:nchar)//'trn', status='replace')
    write(22, '(3E20.12)') xorig, yorig, zorig
    close(22)

    !! Project the surface points to 2-D for triangulation.
    allocate(surface_x(nsurf), surface_y(nsurf), surface_z(nsurf))
    do i = 1, nsurf
       surface_x(i) = all_points(i,1)
       surface_y(i) = all_points(i,2)
       surface_z(i) = all_points(i,3)
    end do

    !! Build a triangulation of the top surface in memory.
    call build_surface_triangulation(nsurf, nedge, edge_seq, surface_x, surface_y, tri_v1, tri_v2, tri_v3, ntri)

    !! Insert non-boundary surface control points into the triangulation so that the
    !! surface mesh still contains every top-surface control point.
    do i = 1, nsurf
       if (tflag(i) == 1) then
          call insert_point_into_triangulation(i, nsurf, surface_x, surface_y, tri_v1, tri_v2, tri_v3, ntri)
       end if
    end do

    !! Write the TetGen PLC.
    !!
    !! The new file is written directly from the geometry in memory.  That removes the
    !! intermediate Triangle output files that previously had to be re-read and matched
    !! back to the cfg points.
    mesh_poly = mesh_cfg
    do i = 1, len_trim(mesh_poly)
       if (mesh_poly(i:i) == '.') then
          mesh_poly(i+1:i+4) = 'poly'
          npre = i
          exit
       end if
    end do
    if (npre == 0) then
       mesh_poly = trim(mesh_cfg)//'.poly'
       npre = len_trim(mesh_cfg) + 1
    end if

    open(11, file=mesh_poly, status='replace', action='write')

    !! Vertex block.
    !!
    !! Top surface vertices first, then the duplicated bottom boundary vertices, then
    !! the interior refinement points.
    write(11, *) nsurf + nedge + nint, 3, 1, 1, ' # Total nodes, dimensions, attributes, boundary marker flag'

    write(11, *) '# Top surface vertices'
    do i = 1, nsurf
       if (tflag(i) == 2) then
          writit = 1
       else
          writit = 0
       end if
       write(11, 113) i, all_points(i,1:3), 1, writit
    end do

    write(11, *) '# Lower boundary vertices'
    do i = 1, nedge
       write(11, 113) nsurf + i, all_points(edge_seq(i),1:2), min_elev, 1, 2
    end do

    if (nint > 0) then
       write(11, *) '# Interior refinement vertices'
       do i = 1, nint
          write(11, 113) nsurf + nedge + i, all_points(nsurf + i,1:3), 1, tflag(nsurf + i)
       end do
    end if

    !! Facet count: top triangles + 2 triangles per side wall + one bottom polygon + PLC fan triangles.
    ntop_facets = ntri
    nside_facets = 2 * nedge
    nbnd_facets = 1
    nplc_facets = 0
    do i = 1, nplc
       nplc_facets = nplc_facets + max(0, plc_npts(i) - 2)
    end do
    nfacets = ntop_facets + nside_facets + nbnd_facets + nplc_facets

    write(11, *)
    write(11, *) nfacets, 1, ' # Total number of facets and boundary marker flag'
    write(11, *) '# 1 = top surface'
    write(11, *) '# 2 = side and bottom boundaries'

    !! Top surface facets.  These triangles are the direct replacement for the Triangle
    !! .ele -> .poly conversion in the original code.
    do i = 1, ntri
       write(11, *) '1 0 1'
       write(11, 122) 3, tri_v1(i), tri_v2(i), tri_v3(i)
    end do

    !! Side walls.  We write them as two triangles per edge so the facets stay planar.
    do i = 1, nedge - 1
       top_id = edge_seq(i)
       bottom_id = nsurf + i
       write(11, *) '1 0 2'
       write(11, 122) 3, top_id, bottom_id, bottom_id + 1
       write(11, *) '1 0 2'
       write(11, 122) 3, top_id, bottom_id + 1, edge_seq(i + 1)
    end do
    top_id = edge_seq(nedge)
    bottom_id = nsurf + nedge
    write(11, *) '1 0 2'
    write(11, 122) 3, top_id, bottom_id, nsurf + 1
    write(11, *) '1 0 2'
    write(11, 122) 3, top_id, nsurf + 1, edge_seq(1)

    !! Bottom boundary.  It is planar, so one polygon facet is sufficient.
    write(11, *) '# Bottom boundary facet'
    write(11, *) '1 0 2'
    write(11, '(I6)', advance='no') nedge
    do i = 1, nedge
       write(11, '(I7)', advance='no') nsurf + i
    end do
    write(11, *)
    write(11, *)

    !! Optional user PLCs.  The input file only lists control points, so here we build
    !! a simple fan triangulation for each PLC.  This keeps the geometry explicit and
    !! avoids the old back-and-forth through Triangle's node/element files.
    if (nplc > 0) then
       write(11, *) '# User defined PLC facets'
       do i = 1, nplc
          do j = 1, plc_npts(i)
             local_id = pmap(plc_pts(i,j))
             if (local_id < 1 .or. local_id > nsurf) then
                stop 'PLC point must lie on the surface set.'
             end if
          end do
          do j = 2, plc_npts(i) - 1
             write(11, *) '1 0 ', plc_bnd(i)
             write(11, 122) 3, pmap(plc_pts(i,1)), pmap(plc_pts(i,j)), pmap(plc_pts(i,j+1))
          end do
       end do
    end if

    write(11, *)
    write(11, *) nholes, ' # Number of holes'
    if (nholes > 0) then
       do i = 1, nholes
          write(11, *) i, holes_xyz(i,:)
       end do
    end if

    write(11, *) nzones, ' # Number of zones'
    if (nzones > 0) then
       do i = 1, nzones
          if (i_flag) then
             write(11, *) zone_labs(i,1), zones_xyz(i,:), zone_labs(i,2), hole_vols(i), zone_sig(i), zone_sigi(i)
          else
             write(11, *) zone_labs(i,1), zones_xyz(i,:), zone_labs(i,2), hole_vols(i)
          end if
       end do
    end if

    close(11)

    !! If TetGen is requested, remove stale output and run the executable.
    inquire(file=mesh_poly(1:npre)//'1.node', exist=exst)
    if (exst) call system('del /Q ' // trim(mesh_poly(1:npre)//'1.node'))
    inquire(file=mesh_poly(1:npre)//'1.ele', exist=exst)
    if (exst) call system('del /Q ' // trim(mesh_poly(1:npre)//'1.ele'))
    inquire(file=mesh_poly(1:npre)//'1.face', exist=exst)
    if (exst) call system('del /Q ' // trim(mesh_poly(1:npre)//'1.face'))
    inquire(file=mesh_poly(1:npre)//'1.neigh', exist=exst)
    if (exst) call system('del /Q ' // trim(mesh_poly(1:npre)//'1.neigh'))

    if (tetflag /= 0) then
       write(*, *) 'Calling TetGen from BuildMesh2'
       if (nzones > 0) then
          call system(trim(tetcom)//' -pnq'//trim(qual)//'a'//trim(max_vol)//'aAA ' // trim(mesh_poly))
       else
          call system(trim(tetcom)//' -pnq'//trim(qual)//'a'//trim(max_vol)//'AA ' // trim(mesh_poly))
       end if
    end if

    write(*, *) 'DONE BUILDING ', trim(mesh_poly)

  contains

    !!----------------------------------------------------------------------------------------------
    subroutine build_surface_triangulation(nsurf_local, nedge_local, edge_seq_local, x, y, tri1, tri2, tri3, ntri_out)
      !! Triangulate the top surface in XY space.
      !!
      !! The boundary loop is triangulated first with a simple ear-clipping pass.
      !! Surface control points with flag=1 are then inserted into the active triangle
      !! list one at a time, which keeps every top-surface control point on the final
      !! surface without relying on Triangle's external output files.

      implicit none

      integer, intent(in) :: nsurf_local, nedge_local
      integer, intent(in) :: edge_seq_local(nedge_local)
      real*8, intent(in) :: x(nsurf_local), y(nsurf_local)
      integer, allocatable, intent(out) :: tri1(:), tri2(:), tri3(:)
      integer, intent(out) :: ntri_out

      integer :: i, ip, inx, nv, nleft, t, pos, j
      integer, allocatable :: poly(:), prev(:), next(:)
      logical, allocatable :: alive(:)
      logical :: ear, ccw
      real*8 :: ax, ay, bx, by, cx, cy, px, py
      real*8 :: cross_val, poly_area
      integer :: maxtri

      !! The top surface cannot exist without a boundary loop.
      if (nedge_local < 3) stop 'Surface triangulation requires at least three boundary points.'

      allocate(poly(nedge_local), prev(nedge_local), next(nedge_local), alive(nedge_local))
      do i = 1, nedge_local
         poly(i) = edge_seq_local(i)
         prev(i) = i - 1
         if (prev(i) < 1) prev(i) = nedge_local
         next(i) = i + 1
         if (next(i) > nedge_local) next(i) = 1
         alive(i) = .true.
      end do

      !! Make sure the boundary loop is counter-clockwise in XY.
      poly_area = 0.0d0
      do i = 1, nedge_local
         j = i + 1
         if (j > nedge_local) j = 1
         poly_area = poly_area + x(poly(i))*y(poly(j)) - x(poly(j))*y(poly(i))
      end do
      ccw = (poly_area > 0.0d0)
      if (.not. ccw) then
         do i = 1, nedge_local/2
            j = poly(i)
            poly(i) = poly(nedge_local - i + 1)
            poly(nedge_local - i + 1) = j
         end do
      end if

      !! Allocate enough storage for the boundary triangulation plus point insertions.
      maxtri = max(3 * nsurf_local, 12)
      allocate(tri1(maxtri), tri2(maxtri), tri3(maxtri))
      ntri_out = 0

      !! Ear clipping on the boundary polygon.
      nleft = nedge_local
      do while (nleft > 3)
         ear = .false.
         do i = 1, nedge_local
            if (.not. alive(i)) cycle
            ip = prev(i)
            inx = next(i)
            if (.not. alive(ip) .or. .not. alive(inx)) cycle

            ax = x(poly(ip)); ay = y(poly(ip))
            bx = x(poly(i));  by = y(poly(i))
            cx = x(poly(inx)); cy = y(poly(inx))
            cross_val = cross2(ax, ay, bx, by, cx, cy)
            if (cross_val <= 0.0d0) cycle

            !! Reject ears that contain any other active boundary vertex.
            ear = .true.
            do j = 1, nedge_local
               if (.not. alive(j)) cycle
               if (j == i .or. j == ip .or. j == inx) cycle
               px = x(poly(j))
               py = y(poly(j))
               if (point_in_triangle(px, py, ax, ay, bx, by, cx, cy)) then
                  ear = .false.
                  exit
               end if
            end do

            if (ear) then
               ntri_out = ntri_out + 1
               tri1(ntri_out) = poly(ip)
               tri2(ntri_out) = poly(i)
               tri3(ntri_out) = poly(inx)
               alive(i) = .false.
               next(ip) = inx
               prev(inx) = ip
               nleft = nleft - 1
               exit
            end if
         end do

         if (.not. ear) stop 'Unable to clip an ear from the boundary polygon.'
      end do

      !! Add the final triangle left in the polygon.
      do i = 1, nedge_local
         if (alive(i)) then
            ntri_out = ntri_out + 1
            tri1(ntri_out) = poly(i)
            tri2(ntri_out) = poly(next(i))
            tri3(ntri_out) = poly(prev(i))
         end if
      end do

    end subroutine build_surface_triangulation

    !!----------------------------------------------------------------------------------------------
    subroutine insert_point_into_triangulation(pt_idx, nsurf_local, x, y, tri1, tri2, tri3, ntri_local)
      !! Insert a single interior surface point into the current triangle list.
      !!
      !! This is not a full Delaunay refinement step.  The goal is more modest: keep
      !! the surface mesh explicit, deterministic, and self-contained, while avoiding
      !! the fragile read-back of Triangle output files.

      implicit none

      integer, intent(in) :: pt_idx, nsurf_local
      real*8, intent(in) :: x(nsurf_local), y(nsurf_local)
      integer, allocatable, intent(inout) :: tri1(:), tri2(:), tri3(:)
      integer, intent(inout) :: ntri_local

      integer :: t, old_ntri, keep_idx, a, b, c
      integer, allocatable :: new1(:), new2(:), new3(:)
      real*8 :: ax, ay, bx, by, cx, cy, px, py

      px = x(pt_idx)
      py = y(pt_idx)
      do t = 1, ntri_local
         a = tri1(t)
         b = tri2(t)
         c = tri3(t)
         ax = x(a); ay = y(a)
         bx = x(b); by = y(b)
         cx = x(c); cy = y(c)
         if (point_in_triangle(px, py, ax, ay, bx, by, cx, cy)) then
            !! Replace one triangle with three triangles that include the new point.
            old_ntri = ntri_local
            allocate(new1(size(tri1) + 2), new2(size(tri2) + 2), new3(size(tri3) + 2))
            if (t > 1) then
               new1(1:t-1) = tri1(1:t-1)
               new2(1:t-1) = tri2(1:t-1)
               new3(1:t-1) = tri3(1:t-1)
            end if
            new1(t) = a; new2(t) = b; new3(t) = pt_idx
            new1(t+1) = b; new2(t+1) = c; new3(t+1) = pt_idx
            new1(t+2) = c; new2(t+2) = a; new3(t+2) = pt_idx
            if (t < old_ntri) then
               new1(t+3:old_ntri+2) = tri1(t+1:old_ntri)
               new2(t+3:old_ntri+2) = tri2(t+1:old_ntri)
               new3(t+3:old_ntri+2) = tri3(t+1:old_ntri)
            end if
            deallocate(tri1, tri2, tri3)
            allocate(tri1(old_ntri + 2), tri2(old_ntri + 2), tri3(old_ntri + 2))
            tri1 = new1(1:old_ntri+2)
            tri2 = new2(1:old_ntri+2)
            tri3 = new3(1:old_ntri+2)
            deallocate(new1, new2, new3)
            ntri_local = old_ntri + 2
            return
         end if
      end do

      !! If we get here, the point could not be placed cleanly.
      !! A warning is safer than silently dropping the point from the surface.
      write(*, *) 'Warning: could not insert surface point ', pt_idx, ' into the triangulation.'

    end subroutine insert_point_into_triangulation

    !!----------------------------------------------------------------------------------------------
    real*8 function cross2(ax, ay, bx, by, cx, cy)
      !! 2-D cross product of the vectors AB and AC.

      implicit none
      real*8, intent(in) :: ax, ay, bx, by, cx, cy
      cross2 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    end function cross2

    !!----------------------------------------------------------------------------------------------
    logical function point_in_triangle(px, py, ax, ay, bx, by, cx, cy)
      !! Barycentric point-in-triangle test in XY space.

      implicit none
      real*8, intent(in) :: px, py, ax, ay, bx, by, cx, cy
      real*8 :: den, s, t, u, tol2

      tol2 = 1.0d-12
      den = cross2(ax, ay, bx, by, cx, cy)
      if (abs(den) < tol2) then
         point_in_triangle = .false.
         return
      end if

      s = cross2(px, py, bx, by, cx, cy) / den
      t = cross2(ax, ay, px, py, cx, cy) / den
      u = 1.0d0 - s - t
      point_in_triangle = (s >= -tol2 .and. t >= -tol2 .and. u >= -tol2)
    end function point_in_triangle

  end subroutine build_tetgen_direct

!__________________________________________________________________________________________________

  113 format(I10,3ES22.12,2I10)
  122 format(4I10)

end module buildmesh2