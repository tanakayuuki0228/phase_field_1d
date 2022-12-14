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

    !________________________________________________________________
    !フェーズフィールド方程式の左辺係数行列を計算するサブルーチン
    subroutine calc_K(K,u)
        double precision,intent(in out),dimension(:) :: K !全体方程式のK
        double precision,intent(in),dimension(:) :: u !変位
        double precision :: K_11,K_12 !要素方程式のK
        !配列初期化
        K=0
        do ele=1,num_ele
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
    !フェーズフィールドの方程式の右辺ベクトルを計算するサブルーチン
    subroutine calc_Bc(Bc)
        double precision,intent(in out),dimension(:) :: Bc !全体の方程式のBc
        double precision :: Bc1 !要素方程式のBc
        !配列の初期化
        Bc=0
        !要素方程式のBcを計算
        Bc1=Gc*dx/2d0/L_0
        !Bcの記憶配列へ格納
        do ele=1,num_ele
            Bc(Bc_num(ele,1))=Bc(Bc_num(ele,1))+Bc1
            Bc(Bc_num(ele,2))=Bc(Bc_num(ele,2))+Bc1
        end do
    end subroutine calc_Bc
    !================================================================
end module matrix