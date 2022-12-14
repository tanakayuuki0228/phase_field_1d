module parameters
    !================================================================
    !このモジュールには
    !・数値的/材料パラメーターの設定
    !・変位規定境界に用いるパラメーターの設定
    !・Gauss関数
    !・CFL条件の確認
    !が含まれる
    !================================================================

    implicit none
    !================================================================
    !材料パラメーター(有次元)
    !================================================================
    !質量密度[kg/m**3]
    double precision,parameter :: density=2450d0
    !縦弾性係数[Pa]
    double precision,parameter :: E=32d0*10d0**9d0
    !ポアソン比
    double precision,parameter :: nyu=0.2d0
    !第一ラメ定数[Pa]
    double precision,parameter :: lamda=nyu*E/(1d0+nyu)/(1d0-2d0*nyu)
    !第二ラメ定数[Pa]
    double precision,parameter :: myu=E/2d0/(1d0+nyu)
    !縦波伝搬速度[m/s]
    double precision,parameter :: v_d=((lamda+2d0*myu)/density)**0.5d0
    !横波伝搬速度[m/s]
    double precision,parameter :: v_s=(myu/density)**0.5d0
    !レーリー波伝搬速度[m/s]
    double precision,parameter :: v_R=(0.862d0+1.14d0*nyu)*v_s/(1d0+nyu)
    !Griffithの破壊エネルギー[J/m**2]
    double precision,parameter :: Gc=3d0
    !================================================================    
    
    
    
    !================================================================
    !変位規定境界に用いるパラメーター
    !================================================================
    !sin波の振幅[m]
    double precision,parameter :: Amp=1.0d0/10d0**9d0
    !sin波の周期[s]
    double precision,parameter :: period=1.0d0/10d0**5d0
    !円周率
    double precision,parameter :: pi=3.1425926589d0
    !角振動数[/s]
    double precision,parameter :: w=2d0*pi/period
    !波長[m]
    double precision,parameter :: wavelen=v_d*period
    !ガウス関数の立ち上がり時間を決めるパラメーター[s]
    double precision,parameter :: Gau=2d0*period
    !================================================================
    
    
    
    !================================================================
    !数値的パラメーター(有次元)
    !================================================================
    !x方向全長[m]
    double precision,parameter :: L_x=0.1d0
    !亀裂近似幅に影響する微小パラメーター[m]
    double precision,parameter :: L_0=0.01d0
    !x方向メッシュサイズ[m]（代表長）
    double precision,parameter :: dx=L_0/5d0
    !時間ステップ幅[s]
    double precision,parameter :: dt=dx/v_d
    !総時間ステップ数
    integer,parameter :: total_timestep=25000
    !解析される総時間[s]
    double precision,parameter :: analyzed_time=dt*total_timestep
    !================================================================
    


    contains
    !================================================================
    !Gauss関数
    !================================================================
    function Gauss(t,n)
        double precision,intent(in) :: t
        integer,intent(in) :: n
        double precision :: Gauss
        double precision :: t_0=3d0*Gau
        if(n==0) then
            !if(t<2d0*t_0) then
                Gauss=exp(-(t-t_0)**2d0/2d0/Gau**2d0)
            !else
                !Gauss=0
            !end if
        else if(n==1) then
            !if(t<2d0*t_0) then
                Gauss=-(t-t_0)/Gau**2d0*exp(-(t-t_0)**2d0/2d0/Gau**2d0)
            !else
                !Gauss=0
            !end if
        else
            !if(t<2d0*t_0) then
                Gauss=-(1d0-(t-t_0)**2d0/Gau**2)/Gau**2d0*exp(-(t-t_0)**2d0/2d0/Gau**2d0)
            !else
                !Gauss=0
            !end if
        end if
    end function Gauss
    !================================================================
    
    
    
    !================================================================
    !CFL条件を確認するサブルーチン
    !================================================================
    subroutine CFL(dx,dt)
        double precision,intent(in) :: dx,dt
        double precision :: cfl_value
        cfl_value=dx/v_d-dt
        if(cfl_value<0) then
            write(*,*) "CFL condition is violated"
        else
            write(*,*) "CFL condition is cleared"
        end if
    end subroutine CFL
    !================================================================
    
end module parameters