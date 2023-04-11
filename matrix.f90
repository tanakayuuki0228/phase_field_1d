module matrix
    !================================================================
    !このモジュールには
    !・方程式の左辺係数行列，右辺ベクトルの計算
    !が含まれる．
    !================================================================

    !================================================================
    !このモジュールで使用される他モジュール
    use parameters !数値的/材料パラメーターを参照
    use mesh !要素数,係数行列/ベクトルの記憶配列,要素と記憶配列紐付ける関数を参照
    !================================================================

    implicit none
    !================================================================
    integer,private :: ele !要素番号ele
    !================================================================

    contains
    !================================================================
    !各方程式の左辺係数行列，右辺ベクトルを計算するサブルーチン
    !================================================================
    !________________________________________________________________
    !lumped質量行列を計算するサブルーチン
    subroutine calc_M(M)
        double precision,intent(in out),dimension(:) :: M
        M(1)=density*dx*5d0/6d0
        M(2)=density*dx
    end subroutine calc_M
    !________________________________________________________________
    !non-lumed質量行列を計算するサブルーチン
    ! subroutine calc_M(M)
    !     double precision,intent(in out),dimension(:) :: M
    !     do ele=1,num_nod-3
    !         M(2*ele-1)=2d0*density*dx/3d0
    !         M(2*ele)=density*dx/6d0
    !     end do
    !     M(2*(num_nod-2)-1)=2d0*density*dx/3d0
    ! end subroutine calc_M
    !________________________________________________________________
    ! !フェーズフィールド方程式の左辺係数行列を計算するサブルーチン
    ! subroutine calc_K(K,u)
    !     double precision,intent(in out),dimension(:) :: K !全体方程式のK
    !     double precision,intent(in),dimension(:) :: u !変位
    !     double precision :: K_11,K_12 !要素方程式のK
    !     !配列初期化
    !     K=0
    !     do ele=1,nint(num_ele/2d0)-10
    !         !要素ごとの(c=1とした時の)ひずみエネルギーを計算
    !         psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
    !         !要素方程式のKの値を計算
    !         K_11=(2d0*psi(ele)+Gc/L_0)*dx/3d0+Gc*L_0/dx
    !         K_12=(2d0*psi(ele)+Gc/L_0)*dx/6d0-Gc*L_0/dx
    !         !Kの記憶配列へ格納
    !         K(K_num(ele,1,1))=K(K_num(ele,1,1))+K_11
    !         K(K_num(ele,1,2))=K(K_num(ele,1,2))+K_12
    !         K(K_num(ele,2,2))=K(K_num(ele,2,2))+K_11
    !     end do
    !     do ele=nint(num_ele/2d0)-9,nint(num_ele/2d0)+10
    !         !要素ごとの(c=1とした時の)ひずみエネルギーを計算
    !         psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
    !         !要素方程式のKの値を計算
    !         K_11=(2d0*psi(ele)+Gc_weak/L_0)*dx/3d0+Gc_weak*L_0/dx
    !         K_12=(2d0*psi(ele)+Gc_weak/L_0)*dx/6d0-Gc_weak*L_0/dx
    !         !Kの記憶配列へ格納
    !         K(K_num(ele,1,1))=K(K_num(ele,1,1))+K_11
    !         K(K_num(ele,1,2))=K(K_num(ele,1,2))+K_12
    !         K(K_num(ele,2,2))=K(K_num(ele,2,2))+K_11
    !     end do
    !     do ele=nint(num_ele/2d0)+11,num_ele
    !         !要素ごとの(c=1とした時の)ひずみエネルギーを計算
    !         psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
    !         !要素方程式のKの値を計算
    !         K_11=(2d0*psi(ele)+Gc/L_0)*dx/3d0+Gc*L_0/dx
    !         K_12=(2d0*psi(ele)+Gc/L_0)*dx/6d0-Gc*L_0/dx
    !         !Kの記憶配列へ格納
    !         K(K_num(ele,1,1))=K(K_num(ele,1,1))+K_11
    !         K(K_num(ele,1,2))=K(K_num(ele,1,2))+K_12
    !         K(K_num(ele,2,2))=K(K_num(ele,2,2))+K_11
    !     end do
    ! end subroutine calc_K
    !________________________________________________________________
    !フェーズフィールド方程式の左辺係数行列を計算するサブルーチン 両端c=1
    subroutine calc_K(K,u)
        double precision,intent(in out),dimension(:) :: K !全体方程式のK
        double precision,intent(in),dimension(:) :: u !変位
        double precision :: K_11,K_12 !要素方程式のK
        !配列初期化
        K=0
        ele=1
        psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
        K_11=(2d0*psi(ele)+Gc/L_0)*dx/3d0+Gc*L_0/dx
        K(K_num(ele,2,2))=K(K_num(ele,2,2))+K_11
        do ele=2,num_ele-1
            !要素ごとの(c=1とした時の)ひずみエネルギーを計算
            psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
            !要素方程式のKの値を計算
            K_11=(2d0*psi(ele)+Gc/L_0)*dx/3d0+Gc*L_0/dx
            K_12=(2d0*psi(ele)+Gc/L_0)*dx/6d0-Gc*L_0/dx
            !Kの記憶配列へ格納
            K(K_num(ele,1,1))=K(K_num(ele,1,1))+K_11
            K(K_num(ele,1,2))=K(K_num(ele,1,2))+K_12
            K(K_num(ele,2,2))=K(K_num(ele,2,2))+K_11
        end do
        ele=num_ele
        psi(ele)=(lamda/2d0+myu)*((u(ele+1)-u(ele))/dx)**2d0
        K_11=(2d0*psi(ele)+Gc/L_0)*dx/3d0+Gc*L_0/dx
        K(K_num(ele,1,1))=K(K_num(ele,1,1))+K_11
    end subroutine calc_K
    !________________________________________________________________

    !________________________________________________________________
    !変位方程式の右辺ベクトルを計算するサブルーチン
    subroutine calc_Bu(Bu,u,a,c)
        double precision,intent(in out),dimension(:) :: Bu !全体の方程式のBu
        double precision,intent(in),dimension(:) :: u,a,c !変位，加速度,フェーズフィールド
        double precision :: Bu_ele !要素方程式のBu
        !配列初期化
        Bu=0
        do ele=1,num_ele
            !一応初期化
            Bu_ele=0
            !要素ごとの応力を計算
            sigma(ele)=(c(ele+1)**2d0+c(ele+1)*c(ele)+c(ele)**2d0)/3d0*(lamda+2d0*myu)*(u(ele+1)-u(ele))/dx
            !要素方程式のBuを計算
            Bu_ele=sigma(ele)
            !Buの記憶配列へ格納
            if(ele==1) then
                Bu(Bu_num(ele,2))=Bu(Bu_num(ele,2))-Bu_ele
            else if(ele<num_ele) then
                Bu(Bu_num(ele,1))=Bu(Bu_num(ele,1))+Bu_ele
                Bu(Bu_num(ele,2))=Bu(Bu_num(ele,2))-Bu_ele
            else !(ele=num_ele)
                Bu(Bu_num(ele,1))=Bu(Bu_num(ele,1))+Bu_ele-density*dx*a(num_nod)/6d0
            end if
        end do
    end subroutine calc_Bu
    !________________________________________________________________

    !________________________________________________________________
    ! !フェーズフィールドの方程式の右辺ベクトルを計算するサブルーチン
    ! subroutine calc_Bc(Bc)
    !     double precision,intent(in out),dimension(:) :: Bc !全体の方程式のBc
    !     double precision :: Bc1 !要素方程式のBc
    !     !配列の初期化
    !     Bc=0
    !     !要素方程式のBcを計算
    !     Bc1=Gc*dx/2d0/L_0
    !     !Bcの記憶配列へ格納
    !     do ele=1,nint(num_ele/2d0)-10
    !         Bc(Bc_num(ele,1))=Bc(Bc_num(ele,1))+Bc1
    !         Bc(Bc_num(ele,2))=Bc(Bc_num(ele,2))+Bc1
    !     end do
    !     do ele=nint(num_ele/2d0)+11,num_ele
    !         Bc(Bc_num(ele,1))=Bc(Bc_num(ele,1))+Bc1
    !         Bc(Bc_num(ele,2))=Bc(Bc_num(ele,2))+Bc1
    !     end do
    !     Bc1=Gc_weak*dx/2d0/L_0
    !     do ele=nint(num_ele/2d0)-9,nint(num_ele/2d0)+10
    !         Bc(Bc_num(ele,1))=Bc(Bc_num(ele,1))+Bc1
    !         Bc(Bc_num(ele,2))=Bc(Bc_num(ele,2))+Bc1
    !     end do
    ! end subroutine calc_Bc
    !________________________________________________________________
    !フェーズフィールドの方程式の右辺ベクトルを計算するサブルーチン 両端c=1
    subroutine calc_Bc(Bc,u)
        double precision,intent(in out),dimension(:) :: Bc !全体の方程式のBc
        double precision,intent(in),dimension(:) :: u !変位
        double precision :: Bc1 !要素方程式のBc
        double precision :: psi_
        double precision :: K_12
        !配列の初期化
        Bc=0
        !要素方程式のBcを計算
        Bc1=Gc*dx/L_0
        !Bcの記憶配列へ格納
        do ele=1,num_nod-2
            Bc(ele)=Bc1
        end do
        psi_=(lamda/2d0+myu)*((u(1+1)-u(1))/dx)**2d0
        K_12=(2d0*psi_+Gc/L_0)*dx/6d0-Gc*L_0/dx
        Bc(1)=Bc(1)-K_12
        psi_=(lamda/2d0+myu)*((u(num_nod)-u(num_nod-1))/dx)**2d0
        K_12=(2d0*psi_+Gc/L_0)*dx/6d0-Gc*L_0/dx
        Bc(num_nod-2)=Bc(num_nod-2)-K_12
    end subroutine calc_Bc
    !================================================================
end module matrix