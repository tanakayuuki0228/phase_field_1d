program main
    !================================================================
    !これはメインプログラムで
    !・PARDISOの設定
    !・変位場/フェーズフィールドの更新計算
    !・計算結果の出力
    !が含まれる．
    !================================================================

    ! include 'C:\Program Files (x86)\Intel\oneAPI\mkl\2022.2.0\include\mkl_pardiso.f90' !intel_mklのPARDISOを使う
    !INCLUDE 'link_fnl_hpc.h'
    !================================================================
    !このモジュール内で使用される他モジュール
    use parameters !時間ステップ幅を参照
    use mesh !節点数を参照
    use variables !変数配列を参照
    use matrix !左辺係数行列,右辺ベクトル計算サブルーチンを参照
    !================================================================

    implicit none
    external pardiso
    !================================================================
    !変数宣言
    integer :: timestep !時間ステップカウンター
    double precision :: t !時刻
    integer :: i !配列番号
    character(57) :: path='C:\Users\tanaka\Documents\phase-field_1d_results\result1\' !出力ファイルのpath
    integer :: record !規定変位境界出力回数
    integer,parameter :: output_interval=100 !何ステップごとに出力するか
    !================================================================
    !pardisoが使う配列
    integer,dimension(64) :: pt
    pt=0 !first call のときに0に初期化されていなくてはならない,途中で自分が書き換えてはいけない
    !================================================================

    !パラメータの出力
    call output_para
    !CFL条件を確認
    call CFL(dx,dt)
    !================================================================
    !計算アルゴリズム
    !================================================================
    !________________________________________________________________
    !変数の初期化
    timestep=0
    t=0
    u=0d0
    v=0d0
    a=0d0
    c=0d0
    record=1
    !________________________________________________________________
    !変位，速度の初期値
    u=0d0
    v=0d0
    !________________________________________________________________
    !変位規定境界の変位，速度，加速度を更新
    call calc_pre_dis
    !________________________________________________________________
    !フェーズフィールドの更新(初期値)
    call allo_K
    call allo_Bc
    call allo_psi
    call calc_phasefield(c,u)
    write(*,*) 'Hello!'
    call output_psi
    call deallo_K
    call deallo_Bc
    call deallo_psi
    !________________________________________________________________
    !加速度の更新(初期値)
    call allo_M
    call allo_Bu
    call allo_sigma
    call calc_acceleration(a,u,c)
    call output_sigma
    call deallo_M
    call deallo_Bu
    call deallo_sigma
    !________________________________________________________________
    !初期値の出力
    call output_u
    call output_v
    call output_a
    call output_c
    call output_predis
    call output_timeshape
    record=2
    !________________________________________________________________
    
    !時間ステップ更新，繰り返し計算
    do timestep=1,total_timestep
        t=dt*timestep
        !____________________________________________________________
        !変位規定境界の変位，速度，加速度を更新
        call calc_pre_dis
        !call output_gauss
        call output_predis
        !____________________________________________________________
        !変位の更新
        call calc_displacement(u,v,a)
        !____________________________________________________________
        !速度の更新(前の時間ステップの加速度の分)
        call calc_velocity(v,a)
        !____________________________________________________________
        !フェーズフィールドの更新
        call allo_K
        call allo_Bc
        call allo_psi
        call calc_phasefield(c,u)
        if(mod(timestep,output_interval)==0) then
        call output_psi
        end if
        call deallo_K
        call deallo_Bc
        call deallo_psi
        !____________________________________________________________
        !加速度の更新
        call allo_M
        call allo_Bu
        call allo_sigma
        call calc_acceleration(a,u,c)
        if(mod(timestep,output_interval)==0) then
            call output_sigma
        end if
        call deallo_M
        call deallo_Bu
        call deallo_sigma
        !____________________________________________________________
        !速度の更新(今の時間ステップの加速度の分)
        call calc_velocity(v,a)
        !____________________________________________________________
        !output_intervalに一度，計算結果の出力
        if(mod(timestep,output_interval)==0) then
            call output_u
            call output_v
            call output_a
            call output_c
            if(mod(timestep,output_interval*100)==0) then
                write(*,*) 'step ',timestep,' is completed'
            end if
        end if
        !____________________________________________________________
        call output_timeshape
    end do
    !================================================================
    

    contains
    !================================================================
    !u,v,a,cを更新するサブルーチン
    !================================================================
    !変位規定境界の変位，速度，加速度を更新するサブルーチン
    subroutine calc_pre_dis
        !(変位規定境界の)変位，速度，加速度を更新
        u(num_nod)=0.5d0*L_x*t**2d0
        v(num_nod)=L_x*t
        a(num_nod)=L_x
    end subroutine calc_pre_dis
    !________________________________________________________________
    !変位を更新するサブルーチン
    subroutine calc_displacement(u,v,a)
        double precision,intent(in),dimension(:) :: v,a
        double precision,intent(in out),dimension(:) :: u
        !陽的ニューマーク法で(変位規定境界以外の)変位を更新
        do i=2,num_nod-1
            u(i)=u(i)+v(i)*dt+a(i)*dt**2d0/2d0
        end do
    end subroutine calc_displacement
    !________________________________________________________________
    !(今か前，どちらかの時間ステップの加速度で)速度を更新するサブルーチン
    subroutine calc_velocity(v,a)
        double precision,intent(in),dimension(:) :: a !今もしくは前の時間ステップの加速度
        double precision,intent(in out),dimension(:) :: v
        !陽的ニューマーク法で(変位規定境界以外の)速度を更新
        do i=2,num_nod-1
            v(i)=v(i)+a(i)*dt/2d0
        end do
    end subroutine
    !________________________________________________________________
    !加速度を更新するサブルーチン
    subroutine calc_acceleration(a,u,c)
        double precision,intent(in),dimension(:) :: u,c
        double precision,intent(in out),dimension(:) :: a
        !lumed質量行列の計算
        call calc_M(M)
        !右辺ベクトルの計算
        call calc_Bu(Bu,u,a,c)
        !(変位規定境界以外の)加速度を更新
        a(2)=Bu(Bu_row(2))/M(1)
        do i=3,num_nod-2
            a(i)=Bu(Bu_row(i))/M(2)
        end do
        a(num_nod-1)=Bu(Bu_row(num_nod-1))/M(1)
    end subroutine calc_acceleration
    !________________________________________________________________
    !フェーズフィールドを更新するサブルーチン
    subroutine calc_phasefield(c,u)
        double precision,intent(in),dimension(:) :: u
        double precision,intent(in out),dimension(:) :: c
        !================================================================
        !PARDISOの設定
        !================================================================
        ! integer,dimension(64),save :: pt !ここでptを宣言するとaccess violation がでる,原因は不明
        integer,save :: maxfct
        integer,save :: mnum
        integer,save :: mtype
        integer,save :: n
        integer,save :: phase
        integer,save,dimension(num_nod+1) :: ia
        integer,save,dimension(2*num_nod-1) :: ja
        integer,save,dimension(num_nod) :: perm
        integer,save :: nrhs
        integer,save,dimension(64) :: iparm
        integer,save :: msglvl
        integer,save :: error
        if(t==0) then
            pt=0
            maxfct=1
            mnum=1
            mtype=-2 !実対称不定値行列,num_nod×num_nod三重対角行列
            n=num_nod
            phase=13
            nrhs=1
            iparm(1)=0 !iparm(1)=0で,iparm(2)-iparm(64)がdefaultに設定される
            msglvl=0

            !ia(i) aのうち,Kのi行で最初の要素のindex
            do i=1,num_nod
                ia(i)=2*i-1
            end do

            !ia(n+1) aの要素の数+1
            ia(num_nod+1)=2*num_nod
            
            !ja(i) a(i)はKの第何列成分か
            do i=1,num_nod
                ja(2*i-1)=i
            end do
            do i=1,num_nod-1
                ja(2*i)=i+1
            end do
        else if(t==dt) then
            phase=23
        else

        end if
        !================================================================
        !左辺係数行列の計算
        call calc_K(K,u)
        !右辺ベクトルの計算
        call calc_Bc(Bc)
        !PARDISOで連立方程式を解く
        call pardiso(pt,maxfct,mnum,mtype,phase,n,K,ia,ja,perm,nrhs,iparm,msglvl,Bc,c,error)
    end subroutine calc_phasefield
    !================================================================
    
    
    !================================================================
    !計算結果をファイルに出力するサブルーチン
    !================================================================
    subroutine output_u
        character(5) :: step
        character(2) :: name='u_'
        character(64) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,u(i)*10d0**9d0
        end do
        close(10)
    end subroutine output_u
    !================================================================
    subroutine output_v
        character(5) :: step
        character(2) :: name='v_'
        character(64) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,v(i)*10d0**9d0
        end do
        close(10)
    end subroutine output_v
    !================================================================
    subroutine output_a
        character(5) :: step
        character(2) :: name='a_'
        character(64) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,a(i)*10d0**9d0
        end do
        close(10)
    end subroutine output_a
    !================================================================
    subroutine output_c
        character(5) :: step
        character(2) :: name='c_'
        character(64) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,c(i)
        end do
        close(10)
    end subroutine output_c
    !================================================================
    subroutine output_Bu
        character(5) :: step
        character(3) :: name='Bu_'
        character(65) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod-2
            x=dx*i
            write(10,'(e24.12,1x,e24.12)') x,Bu(i)
        end do
        close(10)
    end subroutine output_Bu
    !================================================================
    subroutine output_sigma
        character(5) :: step
        character(6) :: name='sigma_'
        character(68) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_ele
            x=dx*(i-0.5d0)
            write(10,'(e24.12,1x,e24.12)') x,sigma(i)
        end do
        close(10)
    end subroutine output_sigma
    !================================================================
    subroutine output_psi
        character(5) :: step
        character(4) :: name='psi_'
        character(66) :: filename
        double precision :: x
        write(step,"(I5.5)") timestep
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_ele
            x=dx*(i-0.5d0)
            write(10,'(e24.12,1x,e24.12)') x,psi(i)
        end do
        close(10)
    end subroutine output_psi
    !================================================================
    subroutine output_gauss
        character(62) :: file_gauss
        character(5) :: g_name='gauss'
        file_gauss=path//g_name
        if(record==1) then
            open(10,file=file_gauss,status='replace')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') t,Gauss(t,0),Gauss(t,1),Gauss(t,2)
            close(10)
        else
            open(10,file=file_gauss,status='old',position='append')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') t,Gauss(t,0),Gauss(t,1),Gauss(t,2)
            close(10)
        end if
    end subroutine output_gauss
    !================================================================
    subroutine output_predis
        character(63) :: file_predis
        character(6) :: predis_name='predis'
        file_predis=path//predis_name
        if(record==1) then
            open(10,file=file_predis,status='replace')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') &
                t*10d0**6d0,u(num_nod)*10d0**9d0,v(num_nod),a(num_nod)
            close(10)
        else
            open(10,file=file_predis,status='old',position='append')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') &
                t*10d0**6d0,u(num_nod)*10d0**9d0,v(num_nod),a(num_nod)
            close(10)
        end if
    end subroutine output_predis
    !================================================================
    subroutine output_timeshape
        character(64) :: file_timeshape
        character(7) :: t1_name='tshape1'
        character(7) :: t2_name='tshape2'
        character(7) :: t3_name='tshape3'
        integer :: tsha1=nint(num_nod/4d0)
        integer :: tsha2=nint(num_nod/2d0)
        integer :: tsha3=nint(3d0*num_nod/4d0)
        
        if(record==1) then
            file_timeshape=path//t1_name
            open(10,file=file_timeshape,status='replace')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha1)*10d0**9d0
            close(10)
            file_timeshape=path//t2_name
            open(10,file=file_timeshape,status='replace')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha2)*10d0**9d0
            close(10)
            file_timeshape=path//t3_name
            open(10,file=file_timeshape,status='replace')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha3)*10d0**9d0
            close(10)
        else
            file_timeshape=path//t1_name
            open(10,file=file_timeshape,status='old',position='append')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha1)*10d0**9d0
            close(10)
            file_timeshape=path//t2_name
            open(10,file=file_timeshape,status='old',position='append')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha2)*10d0**9d0
            close(10)
            file_timeshape=path//t3_name
            open(10,file=file_timeshape,status='old',position='append')
            write(10,'(e24.12,1x,e24.12)') t*10d0**6d0,u(tsha3)*10d0**9d0
            close(10)
        end if
    end subroutine output_timeshape
    !================================================================
    subroutine output_M
        character(58) :: file_M
        character(1) :: M_name='M'
        file_M=path//M_name
        open(10,file=file_M,status='replace')
        write(10,'(e24.12,1x,e24.12)') M(1),M(2)
        close(10)
    end subroutine output_M
    !================================================================
    
    
    !================================================================
    !パラメーターをファイルに出力するサブルーチン
    !================================================================
    subroutine output_para
        character(10) :: name_p='parameters'
        character(67) :: filename_p
        filename_p=path//name_p
        open(10,file=filename_p,status='replace')
        write(10,'(a,e24.12,a)') 'density = ',density,' [kg/m**3]'
        write(10,'(a,e24.12,a)') 'E = ',E,' [Pa]'
        write(10,'(a,f)') 'nyu = ',nyu
        write(10,'(a,e24.12,a)') 'lamda = ',lamda,' [Pa]'
        write(10,'(a,e24.12,a)') 'myu = ',myu,' [Pa]'
        write(10,'(a,e24.12,a)') 'v_d = ',v_d,' [m/s]'
        write(10,'(a,e24.12,a)') 'v_s = ',v_s,' [m/s]'
        write(10,'(a,e24.12,a)') 'v_R = ',v_R,' [m/s]'
        write(10,'(a,e24.12,a)') 'Gc = ',Gc,' [J/m**2]'
        write(10,'(a,e24.12,a)') 'Amp = ',Amp,' [m]'
        write(10,'(a,e24.12,a)') 'period = ',period,' [s]'
        write(10,'(a,e24.12,a)') 'angular_frecuency = ',w,' [/s]'
        write(10,'(a,e24.12,a)') 'wavelength = ',wavelen,' [m]'
        write(10,'(a,e24.12,a)') 'Gau = ',Gau,' [s]'
        write(10,'(a,e24.12,a)') 'L_x = ',L_x,' [m]'
        write(10,'(a,e24.12,a)') 'dx = ',dx,' [m]'
        write(10,'(a,e24.12,a)') 'L_0 = ',L_0,' [m]'
        write(10,'(a,e24.12,a)') 'dt = ',dt,' [s]'
        write(10,'(a,i)') 'total_timestep = ',total_timestep
        write(10,'(a,e24.12,a)') 'total_time = ',analyzed_time,' [s]'
        write(10,'(a,i)') 'number_node = ',num_nod
        write(10,'(a,e24.12)') 'x1 = ',dx*(nint(num_nod/4d0)-1)
        write(10,'(a,e24.12)') 'x2 = ',dx*(nint(num_nod/2d0)-1)
        write(10,'(a,e24.12)') 'x3 = ',dx*(nint(3d0*num_nod/4d0)-1)
        close(10)
    end subroutine
    !================================================================ 
end program main