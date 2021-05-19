module qsort_mod
 
implicit none
 
type group
    integer :: order    ! original order of unsorted data
    real :: value       ! values to be sorted by
		integer :: weight
		integer :: ival
end type group
 
contains
 
recursive subroutine QSort(a,na)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
type (group), dimension(nA), intent(in out) :: A
 
! LOCAL VARIABLES
integer :: left, right
real :: random
real :: pivot
type (group) :: temp
integer :: marker
 
    if (nA > 1) then
 
        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QSort(A(:marker-1),marker-1)
        call QSort(A(marker:),nA-marker+1)
 
    end if
 
end subroutine QSort
 
end module qsort_mod

PROGRAM SOLVER
	USE QSORT_MOD
	IMPLICIT NONE
	
	integer,parameter :: k1 = selected_int_kind(1)
	INTEGER :: ITEM_COUNT,CAPACITY,IVAL,IWEIGHT
	INTEGER :: I,W,IT,IT2,ILEVEL
	INTEGER(8) :: ITEM_COUNT8,CAPACITY8,I8,IVALMIN,IWEIGHTMIN,EVAL,NLEVEL,MAXPROFIT,IB
	TYPE (GROUP),DIMENSION(:),ALLOCATABLE :: ITEMS_SORTED
	INTEGER(kind=k1),DIMENSION(:),ALLOCATABLE :: TAKEN,TAKEN_TEMP,TAKEN_OPT
	LOGICAL,DIMENSION(:),ALLOCATABLE :: YNOT
	INTEGER,DIMENSION(:),ALLOCATABLE :: VALUE,WEIGHT,ENUMERATION,K_TEMP
	INTEGER,DIMENSION(:,:),ALLOCATABLE :: K
	CHARACTER(12) :: FORMATER
	REAL(8) :: AVGWEIGHT
	
	
	
	
	
	!PRINT *, "Hello World"
	
	OPEN(2223,FILE="tmp.data")
	READ(2223,*) ITEM_COUNT,CAPACITY
	ITEM_COUNT8 = ITEM_COUNT
	CAPACITY8 = CAPACITY
	!PRINT*, ITEM_COUNT,CAPACITY
	!PRINT*, ITEM_COUNT8,CAPACITY8
	I = ITEM_COUNT*CAPACITY
	I8 = ITEM_COUNT8*CAPACITY8
	
	
	ALLOCATE(TAKEN(1:ITEM_COUNT8),VALUE(1:ITEM_COUNT8),WEIGHT(1:ITEM_COUNT8))
	ALLOCATE(ENUMERATION(1:ITEM_COUNT8),ITEMS_SORTED(1:ITEM_COUNT8))
	TAKEN = 0
	
	!GOTO 999
	
	DO I = 1,ITEM_COUNT
		ENUMERATION(I) = I
		READ(2223,*) IVAL, IWEIGHT
		WEIGHT(I) = IWEIGHT
		VALUE(I) = IVAL
	 	ITEMS_SORTED(I)%ORDER = I
	 	ITEMS_SORTED(I)%IVAL = IVAL
	 	ITEMS_SORTED(I)%WEIGHT = IWEIGHT
	 	ITEMS_SORTED(I)%VALUE = REAL(IVAL)/REAL(IWEIGHT)
	ENDDO
	CALL QSort(ITEMS_SORTED,ITEM_COUNT)
	CLOSE(2223)
	IVAL = 0
	IVALMIN = MINVAL(VALUE)
	IWEIGHTMIN = MINVAL(WEIGHT)
	AVGWEIGHT = REAL(SUM(WEIGHT)/ITEM_COUNT8)
	WRITE(2234,*) 'IVALMIN',IVALMIN
	WRITE(2234,*) 'IWEIGHTMIN',IWEIGHTMIN	
	WRITE(2234,*) 'AVGWEIGHT',AVGWEIGHT
  !WHERE(WEIGHT(:) < 0.2*AVGWEIGHT) YNOT = .TRUE.

	
	FLUSH(2234)
	IF(CAPACITY8 < 1000000) THEN
	ALLOCATE(K(0:ITEM_COUNT8,IWEIGHTMIN-1:CAPACITY8))
	K = 0
	K(0,:) = 0
	K(:,IWEIGHTMIN-1) = 0
	DO I = 1,ITEM_COUNT8
   K(I,:) = K(I-1,:)
	 DO	W = IWEIGHTMIN,CAPACITY8
		 IF (WEIGHT(I) < W) THEN
			 IB = MAX(W-WEIGHT(I),IWEIGHTMIN-1)
			 K(I,W) = MAX(VALUE(I)+K(I-1,IB),K(I-1,W))
		 END IF
	 ENDDO
	ENDDO
	DEALLOCATE(VALUE)
	!write out K
	!DO W = IWEIGHTMIN-1,CAPACITY8
	!	WRITE(2235,'(8000I8)') (K(I,W),I=0,ITEM_COUNT8)
	!ENDDO
	
	!end write out K
	 
	 I = ITEM_COUNT8
	 W = CAPACITY8
	 
	 Taker: DO
		IF(K(I,W) > K(I-1,W)) THEN
		 TAKEN(I) = 1
		 W = W-WEIGHT(I)
		 IF(W < IWEIGHTMIN) EXIT Taker
		 I = I-1
		 IF(I < 0) EXIT Taker
	 ELSE
		 I = I-1
		 IF(I < 0) EXIT Taker
	 END IF
		 
	 ENDDO Taker
	IVAL = K(ITEM_COUNT8,CAPACITY8)
	 DEALLOCATE(K,WEIGHT)
	 
 ELSE
 	ALLOCATE(K_TEMP(IWEIGHTMIN-1:CAPACITY8))

 	DO I = 1,ITEM_COUNT8
		
	 IF(I == 1) THEN
	 	K_TEMP(:) = 0
	 ELSE
	 	K_TEMP(:) = K(I-1,:)
	 END IF
	 IF(ALLOCATED(K)) DEALLOCATE(K)
 	 ALLOCATE(K(I-1:I,IWEIGHTMIN-1:CAPACITY8))
	 K(:,:) = 0
 	 K(I-1,:) = K_TEMP(:)
	 
 	 DO	W = IWEIGHTMIN,CAPACITY8
 		 IF (WEIGHT(I) < W) THEN
 			 IB = MAX(W-WEIGHT(I),IWEIGHTMIN-1)
 			 K(I,W) = MAX(VALUE(I)+K(I-1,IB),K(I-1,W))
 		 END IF
 	 ENDDO
 	ENDDO
	
	!DO I = IWEIGHTMIN-1,CAPACITY8
	!WRITE(2234,*) K(:,I)		
	!ENDDO

	WRITE(2234,*) 'maxval ',K(ITEM_COUNT8,CAPACITY8)
	 
	 
 	IVAL = 0
 	IWEIGHT = 0
 	DO I = ITEM_COUNT,1,-1
	!	IF(ITEMS_SORTED(I)%WEIGHT > 150000) CYCLE
 		IF(IWEIGHT+ITEMS_SORTED(I)%WEIGHT < CAPACITY8) THEN
			WRITE(2234,'(I9,F10.5,I9,I9)') ITEMS_SORTED(I)%ORDER,ITEMS_SORTED(I)%VALUE,ITEMS_SORTED(I)%IVAL,ITEMS_SORTED(I)%WEIGHT
 			IWEIGHT = IWEIGHT + ITEMS_SORTED(I)%WEIGHT
 			IVAL = IVAL+ ITEMS_SORTED(I)%IVAL
			TAKEN(ITEMS_SORTED(I)%ORDER) = 1
 		END IF
 	ENDDO 
	 END IF
	 

	 

!!	WRITE(*,*) IVAL,IWEIGHT
	
	!official output section

	WRITE(*,*) IVAL,0
	
	FORMATER = ''
	WRITE(FORMATER,*) ITEM_COUNT
	FORMATER = '('//ADJUSTL(FORMATER(2:12))
	FORMATER = TRIM(FORMATER)
	FORMATER  = TRIM(ADJUSTL(FORMATER))//'I6)'
	WRITE(*,FORMATER) TAKEN
	
	DEALLOCATE(TAKEN)
	
	
	999 CONTINUE
	
	
END PROGRAM SOLVER