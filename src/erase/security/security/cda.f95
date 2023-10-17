PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  IMPLICIT NONE
  TYPE (cube)     :: abcube , acube, bcube, refcube, drho
  !
  CALL cube_get (abcube, "../../test/CuCO+/ab.cube")
  CALL cube_get (acube, "../../test/CuCO+/a.cube")
  CALL cube_get (bcube, "../../test/CuCO+/b.cube")
  !
  CALL cube_write (abcube,"ab_test.cube")
  CALL cube_write (acube, "a_test.cube")
  CALL cube_write (bcube, "b_test.cube")
  !
  !print *, 'a cube integration'
  !write(*,'(E11.4)') cube_int(acube)
  !print *, 'b cube integration'
  !write(*,'(E11.4)') cube_int(bcube)
  !print *, 'ab cube integration'
  !write(*,'(E11.4)') cube_int(abcube)
  !
  !refcube = acube + bcube
  !drho    = abcube - refcube ! drho = abcube - refcube
  refcube = cube_add(acube,bcube)
  drho    = cube_sub(abcube,refcube)
  !
  !CALL cube_write (refcube, "refcube.cube")
  !CALL cube_write (drho, "drhocube.cube")
  !
  print *, 'integral of drho', cube_int(drho)
  CALL cube_cdz (drho, "final_result.csv")  



END PROGRAM cda
