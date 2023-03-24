module mesh
    !================================================================
    !このモジュール内には
    !・要素数，節点数の確定
    !・左辺係数行列，右辺ベクトルの動的配列宣言
    !・配列の動的割り付け
    !が含まれる．
    !================================================================

    !================================================================
    !このモジュール内で使用される他モジュール
    use parameters !数値的パラメーターを参照
    !================================================================

    implicit none

    !================================================================
    !要素数，節点数
    !================================================================
    !x方向の要素数，節点数
    integer,parameter :: num_ele=nint(L_x/dx)
    integer,parameter :: num_nod=num_ele+1
    !________________________________________________________________
    !変位規定境界にある節点の数
    integer,parameter,private :: num_prenod=2
    !================================================================


    !================================================================
    !左辺係数行列，右辺ベクトルの記憶配列
    !================================================================
    !________________________________________________________________
    !Mの記憶配列を作成
    !________________________________________________________________
    !Mに要る記憶容量
    integer,parameter :: size_M=2
    !Mの数値を記憶する配列
    double precision,allocatable,dimension(:) :: M
    !________________________________________________________________


    !________________________________________________________________
    !Kの記憶配列を作成
    !________________________________________________________________
    !Kに要る記憶容量
    ! integer,parameter :: size_K=2*num_nod-1
    integer,parameter :: size_K=2*num_nod-5 !両端c=0
    !Kの数値を記憶する配列
    double precision,allocatable,dimension(:) :: K
    !________________________________________________________________

    !________________________________________________________________
    !Buの記憶配列の作成
    !______________________________
    !Buに要る記憶容量
    integer,parameter :: size_Bu=num_nod-num_prenod
    !Buの数値を記憶する配列
    double precision,allocatable,dimension(:) :: Bu
    !________________________________________________________________

    !________________________________________________________________
    !Bcの記憶配列
    !________________________________________________________________
    !Bcに要る記憶容量
    ! integer,parameter :: size_Bc=num_nod
    integer,parameter :: size_Bc=num_nod-2 !両端c=0
    !Bcの数値を記憶する配列
    double precision,allocatable,dimension(:) :: Bc
    !================================================================
    
    !================================================================
    !要素ごとの計算値の記憶配列
    !================================================================
    integer,parameter :: size_sigma=num_ele
    double precision,allocatable,dimension(:) :: sigma
    !________________________________________________________________
    integer,parameter :: size_psi=num_ele
    double precision,allocatable,dimension(:) :: psi

    contains
    !================================================================
    !配列を割り付け,開放するサブルーチン
    !================================================================
    subroutine allo_M
        allocate(M(size_M))
    end subroutine allo_M

    subroutine deallo_M
        deallocate(M)
    end subroutine deallo_M

    subroutine allo_K
        allocate(K(size_K))
    end subroutine allo_K

    subroutine deallo_K
        deallocate(K)
    end subroutine deallo_K

    subroutine allo_Bu
        allocate(Bu(size_Bu))
    end subroutine allo_Bu

    subroutine deallo_Bu
        deallocate(Bu)
    end subroutine deallo_Bu

    subroutine allo_Bc
        allocate(Bc(size_Bc))
    end subroutine allo_Bc

    subroutine deallo_Bc
        deallocate(Bc)
    end subroutine deallo_Bc

    subroutine allo_sigma
        allocate(sigma(size_sigma))
    end subroutine allo_sigma

    subroutine deallo_sigma
        deallocate(sigma)
    end subroutine deallo_sigma

    subroutine allo_psi
        allocate(psi(size_psi))
    end subroutine allo_psi

    subroutine deallo_psi
        deallocate(psi)
    end subroutine deallo_psi
    !================================================================


    !================================================================
    !要素番号，要素方程式の行列の行/列番号と記憶配列の配列番号を結ぶ関数+α
    !================================================================
    !________________________________________________________________
    !K用
    function K_num(ele,row,col)
        !要素番号ele,要素方程式での行/列番号row,col
        integer,intent(in) :: ele,row,col
        !全体の方程式の記憶配列での配列番号
        integer :: K_num
        ! K_num=2*ele-1+row-1+col-1
        K_num=2*ele-1+row-1+col-1-2 !両端c=0
    end function K_num
    !________________________________________________________________
    !Bu用
    !要素と記憶配列紐付ける関数
    function Bu_num(ele,row)
        !要素番号ele,要素方程式での行番号row
        integer,intent(in) :: ele,row
        !全体の方程式の記憶配列での配列番号
        integer :: Bu_num
        Bu_num=ele-2+row
    end function Bu_num

    !加速度配列番号とBu記憶配列番号を合わせる関数
    function Bu_row(row)
        !加速度配列番号
        integer,intent(in) :: row
        !Bu記憶配列番号
        integer :: Bu_row
        Bu_row=row-1
    end function Bu_row
    !________________________________________________________________
    !Bc用
    function Bc_num(ele,row)
        !要素番号ele,要素方程式での行番号row
        integer,intent(in) :: ele,row
        !全体の方程式の記憶配列での配列番号
        integer :: Bc_num
        ! Bc_num=ele-1+row
        Bc_num=ele-1+row-1 !両端c=0
    end function Bc_num
    !================================================================
    
end module