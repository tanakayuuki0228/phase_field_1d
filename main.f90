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
    use,intrinsic :: iso_fortran_env !現在時刻取得に使用（intrinsicが必要かどうかは不明）
    !================================================================

    implicit none
    external pardiso
    !================================================================
    !変数宣言
    integer :: timestep !時間ステップカウンター
    double precision :: t !時刻
    integer :: i !配列番号
    integer,parameter :: int_path_=49 !pathの文字数
    integer,parameter :: int_path=65 !pathの文字数
    character(int_path_) :: path_='C:\Users\tanaka\Documents\phase-field_1d_results\' !出力ファイルのpath
    character(int_path) :: path
    character(256) :: command
    character(16) :: folder
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    integer :: timevalue(8)
    integer :: record !規定変位境界出力回数
    integer :: output_interval !何ステップごとに出力するか
    integer :: display_interval !step数を画面出力する間隔
    double precision :: energy_initial !t=0における運動エネルギー,ひずみエネルギー,破壊エネルギーの和
    double precision :: energy_max !グラフの縦軸用,エネルギーの最大値
    integer :: count_3,count_4,count_5 !c=0となったtimestepを記憶 
    double precision :: u_break_3,u_break_4,u_break_5 !c=0となったときの右端変位
    double precision :: t_break_3,t_break_4,t_break_5 !c=0となったときの時刻
    !================================================================
    !pardisoが使う配列
    integer,dimension(64) :: pt
    pt=0 !first call のときに0に初期化されていなくてはならない,途中で自分が書き換えてはいけない
    !================================================================

    t=0
    !================================================================
    !出力間隔を決定
    do while (10d0**t<total_timestep)
        t=t+1d0
    end do
    if(t>5d0) then
        output_interval=10**nint(t-5d0)
    else
        output_interval=1
    end if
    !================================================================

    !================================================================
    !現在時刻を取得
    call date_and_time(date,time,zone,timevalue)
    !計算結果を入れるfolderを生成
    write(folder,"(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2)") &
    timevalue(1),'.',timevalue(2),'.',timevalue(3),'.',timevalue(5),'.',timevalue(6)
    path=path_//folder
    write(command,*) 'mkdir ',path
    call system(command)
    !================================================================

    !パラメータの出力
    call output_para
    !CFL条件を確認
    call CFL(dx,dt)
    call output_a_exa
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
    c=1d0
    record=1
    display_interval=1
    count_3=0
    count_4=0
    count_5=0
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
    call deallo_K
    call deallo_Bc
    call deallo_psi
    !________________________________________________________________
    !加速度の更新(初期値)
    call allo_M
    call allo_Bu
    call allo_sigma
    call calc_acceleration(a,u,c)
    call deallo_M
    call deallo_Bu
    call deallo_sigma
    !________________________________________________________________
    !初期値の出力
    call output_u
    call output_v
    call output_a
    call output_c
    call output_sigma
    call output_psi_ela
    call output_psi_fra
    call output_psi_kin
    call output_energy
    record=2
    !________________________________________________________________
    
    !時間ステップ更新，繰り返し計算
    do timestep=1,total_timestep
        t=dt*timestep
        !____________________________________________________________
        !変位規定境界の変位，速度，加速度を更新
        call calc_pre_dis
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
        call deallo_K
        call deallo_Bc
        call deallo_psi
        !____________________________________________________________
        !加速度の更新
        call allo_M
        call allo_Bu
        call allo_sigma
        call calc_acceleration(a,u,c)
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
            call output_sigma
            call output_psi_ela
            call output_psi_fra
            call output_psi_kin
            call output_energy
            call break_check
            
            if(display_interval<100) then
                if(t>analyzed_time*display_interval/100d0) then
                    write(*,'(i2,a)') display_interval,'% is completed'
                    display_interval=display_interval+1
                end if
            end if
        end if
        !____________________________________________________________
    end do

    !matlabのスプリクトを作成
    call output_matlab

    !folderにコードをディレクトリごとコピー
    write(command,*) 'cp -r C:\Users\tanaka\Documents\phase_field_1d_newest_code\phase_field_1d ',path
    call system(command)
    !================================================================

    contains
    !================================================================
    !u,v,a,cを更新するサブルーチン
    !================================================================
    !変位規定境界の変位，速度，加速度を更新するサブルーチン
    subroutine calc_pre_dis
        !(変位規定境界の)変位，速度，加速度を更新
        ! if(t<time_stand) then
        !     u(num_nod)=(t-time_stand/2d0/pi*sin(2d0*pi*t/time_stand))/2d0*v_end
        !     v(num_nod)=(1d0-cos(2d0*pi*t/time_stand))/2d0*v_end
        !     a(num_nod)=pi/time_stand*sin(2d0*pi*t/time_stand)*v_end
        ! else
        !     u(num_nod)=(t-time_stand/2d0)*v_end
        !     v(num_nod)=v_end
        !     a(num_nod)=0d0
        ! end if
        u(num_nod)=0.5d0*L_x*t**2d0 !Amp*Gauss(t,0)*sin(w*t)
        v(num_nod)=L_x*t !Amp*(Gauss(t,1)*sin(w*t)+Gauss(t,0)*w*cos(w*t))
        a(num_nod)=L_x !Amp*(Gauss(t,2)*sin(w*t)+2*Gauss(t,1)*w*cos(w*t)-Gauss(t,0)*w**2d0*sin(w*t))
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
    !加速度を更新するサブルーチン(lumped)
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
    !加速度を更新するサブルーチン(non-lumped)
    ! subroutine calc_acceleration(a,u,c)
    !     double precision,intent(in),dimension(:) :: u,c
    !     double precision,intent(inout),dimension(:) :: a
    !     double precision,dimension(num_nod-2) :: a_sol
    !     !PARDISOの設定
    !     integer,save :: maxfct
    !     integer,save :: mnum
    !     integer,save :: mtype
    !     integer,save :: n
    !     integer,save :: phase
    !     integer,save,dimension(num_nod+1) :: ia
    !     integer,save,dimension(2*num_nod-1) :: ja
    !     integer,save,dimension(num_nod) :: perm
    !     integer,save :: nrhs
    !     integer,save,dimension(64) :: iparm
    !     integer,save :: msglvl
    !     integer,save :: error
    !     if(t==0) then
    !         pt=0
    !         maxfct=1
    !         mnum=1
    !         mtype=-2 !実対称不定値行列,(num_nod-2)×(num_nod-2)三重対角行列
    !         n=num_nod-2
    !         phase=13
    !         nrhs=1
    !         iparm(1)=0 !iparm(1)=0で,iparm(2)-iparm(64)がdefaultに設定される
    !         msglvl=0

    !         !ia(i) aのうち,Mのi行で最初の要素のindex
    !         do i=1,num_nod-2
    !             ia(i)=2*i-1
    !         end do

    !         !ia(n+1) aの要素の数+1
    !         ia(num_nod-1)=2*(num_nod-2)
            
    !         !ja(i) a(i)はMの第何列成分か
    !         do i=1,num_nod-2
    !             ja(2*i-1)=i
    !         end do
    !         do i=1,num_nod-3
    !             ja(2*i)=i+1
    !         end do
    !     else if(t==dt) then
    !         phase=23
    !     else

    !     end if
    !     !================================================================
    !     call calc_M(M)
    !     call calc_Bu(Bu,u,a,c)
    !     call pardiso(pt,maxfct,mnum,mtype,phase,n,M,ia,ja,perm,nrhs,iparm,msglvl,Bu,a_sol,error)
    !     do i=2,num_nod-1
    !         a(i)=a_sol(i-1)
    !     end do
    ! end subroutine calc_acceleration
    !________________________________________________________________
    ! !フェーズフィールドを更新するサブルーチン
    ! subroutine calc_phasefield(c,u)
    !     double precision,intent(in),dimension(:) :: u
    !     double precision,intent(in out),dimension(:) :: c
    !     !================================================================
    !     !PARDISOの設定
    !     !================================================================
    !     !! integer,dimension(64),save :: pt !ここでptを宣言するとaccess violation がでる,原因は不明
    !     integer,save :: maxfct
    !     integer,save :: mnum
    !     integer,save :: mtype
    !     integer,save :: n
    !     integer,save :: phase
    !     integer,save,dimension(num_nod+1) :: ia
    !     integer,save,dimension(2*num_nod-1) :: ja
    !     integer,save,dimension(num_nod) :: perm
    !     integer,save :: nrhs
    !     integer,save,dimension(64) :: iparm
    !     integer,save :: msglvl
    !     integer,save :: error
    !     if(t==0) then
    !         pt=0
    !         maxfct=1
    !         mnum=1
    !         mtype=-2 !実対称不定値行列,num_nod×num_nod三重対角行列
    !         n=num_nod
    !         phase=13
    !         nrhs=1
    !         iparm(1)=0 !iparm(1)=0で,iparm(2)-iparm(64)がdefaultに設定される
    !         msglvl=0

    !         !ia(i) aのうち,Kのi行で最初の要素のindex
    !         do i=1,num_nod
    !             ia(i)=2*i-1
    !         end do

    !         !ia(n+1) aの要素の数+1
    !         ia(num_nod+1)=2*num_nod
            
    !         !ja(i) a(i)はKの第何列成分か
    !         do i=1,num_nod
    !             ja(2*i-1)=i
    !         end do
    !         do i=1,num_nod-1
    !             ja(2*i)=i+1
    !         end do
    !     else if(t==dt) then
    !         phase=23
    !     else

    !     end if
    !     !================================================================
    !     !左辺係数行列の計算
    !     call calc_K(K,u)
    !     !右辺ベクトルの計算
    !     call calc_Bc(Bc)
    !     !PARDISOで連立方程式を解く
    !     call pardiso(pt,maxfct,mnum,mtype,phase,n,K,ia,ja,perm,nrhs,iparm,msglvl,Bc,c,error)
    ! end subroutine calc_phasefield
    !________________________________________________________________
    !フェーズフィールドを更新するサブルーチン 両端c=1
    subroutine calc_phasefield(c,u)
        double precision,intent(in),dimension(:) :: u
        double precision,intent(in out),dimension(:) :: c
        double precision,dimension(num_nod-2) :: c_sol
        !================================================================
        !PARDISOの設定
        !================================================================
        !! integer,dimension(64),save :: pt !ここでptを宣言するとaccess violation がでる,原因は不明
        integer,save :: maxfct
        integer,save :: mnum
        integer,save :: mtype
        integer,save :: n
        integer,save :: phase
        integer,save,dimension(num_nod-1) :: ia
        integer,save,dimension(2*num_nod-5) :: ja
        integer,save,dimension(num_nod-2) :: perm
        integer,save :: nrhs
        integer,save,dimension(64) :: iparm
        integer,save :: msglvl
        integer,save :: error
        if(t==0) then
            pt=0
            maxfct=1
            mnum=1
            mtype=-2 !実対称不定値行列,num_nod×num_nod三重対角行列
            n=num_nod-2
            phase=13
            nrhs=1
            iparm(1)=0 !iparm(1)=0で,iparm(2)-iparm(64)がdefaultに設定される
            msglvl=0

            !ia(i) aのうち,Kのi行で最初の要素のindex
            do i=1,num_nod-2
                ia(i)=2*i-1
            end do

            !ia(n+1) aの要素の数+1
            ia(num_nod-1)=2*num_nod-4
            
            !ja(i) a(i)はKの第何列成分か
            do i=1,num_nod-2
                ja(2*i-1)=i
            end do
            do i=1,num_nod-3
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
        call calc_Bc(Bc,u)
        !PARDISOで連立方程式を解く
        call pardiso(pt,maxfct,mnum,mtype,phase,n,K,ia,ja,perm,nrhs,iparm,msglvl,Bc,c_sol,error)
        do i=2,num_nod-1
            c(i)=c_sol(i-1)
        end do
    end subroutine calc_phasefield
    !================================================================
    
    
    !================================================================
    !計算結果をファイルに出力するサブルーチン
    !================================================================
    subroutine output_u
        character(5) :: step
        character(5) :: name='\u\u_'
        character(int_path+10) :: filename
        double precision :: x
        write(command,*) 'mkdir ',path,'\u'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,u(i)
        end do
        close(10)
    end subroutine output_u
    !================================================================
    subroutine output_v
        character(5) :: step
        character(5) :: name='\v\v_'
        character(int_path+10) :: filename
        double precision :: x
        write(command,*) 'mkdir ',path,'\v'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,v(i)
        end do
        close(10)
    end subroutine output_v
    !================================================================
    subroutine output_a
        character(5) :: step
        character(5) :: name='\a\a_'
        character(int_path+10) :: filename
        double precision :: x
        write(command,*) 'mkdir ',path,'\a'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,a(i)
        end do
        close(10)
    end subroutine output_a
    !================================================================
    subroutine output_c
        character(5) :: step
        character(5) :: name='\c\c_'
        character(int_path+10) :: filename
        double precision :: x
        write(command,*) 'mkdir ',path,'\c'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            write(10,'(e24.12,1x,e24.12e4)') x,c(i)
        end do
        close(10)
    end subroutine output_c
    !================================================================
    subroutine output_sigma
        character(5) :: step
        character(9) :: name='\sig\sig_'
        character(int_path+14) :: filename
        double precision :: x,sig
        write(command,*) 'mkdir ',path,'\sig'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_ele
            x=dx*(i-1)
            sig=c(i)**2d0*(lamda+2d0*myu)*(u(i+1)-u(i))/dx
            write(10,'(e24.12,1x,e24.12e4)') x,sig
            x=dx*i
            sig=c(i+1)**2d0*(lamda+2d0*myu)*(u(i+1)-u(i))/dx
            write(10,'(e24.12,1x,e24.12e4)') x,sig
        end do
        close(10)
    end subroutine output_sigma
    !================================================================
    subroutine output_psi_ela
        character(5) :: step
        character(16) :: name='\psi_ela\psi_ela_'
        character(int_path+21) :: filename
        double precision :: x,psi_ela
        write(command,*) 'mkdir ',path,'\psi_ela'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_ele
            x=dx*(i-1)
            psi_ela=c(i)**2d0*(lamda/2d0+myu)*((u(i+1)-u(i))/dx)**2d0
            write(10,'(e24.12,1x,e24.12e4)') x,psi_ela
            x=dx*i
            psi_ela=c(i+1)**2d0*(lamda/2d0+myu)*((u(i+1)-u(i))/dx)**2d0
            write(10,'(e24.12,1x,e24.12e4)') x,psi_ela
        end do
        close(10)
    end subroutine output_psi_ela
    !================================================================
    subroutine output_psi_kin
        character(5) :: step
        character(16) :: name='\psi_kin\psi_kin_'
        character(int_path+21) :: filename
        double precision :: x,psi_kin
        write(command,*) 'mkdir ',path,'\psi_kin'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_nod
            x=dx*(i-1)
            psi_kin=density*v(i)**2d0/2d0
            write(10,'(e24.12,1x,e24.12e4)') x,psi_kin
        end do
        close(10)
    end subroutine output_psi_kin
    !================================================================
    subroutine output_psi_fra
        character(5) :: step
        character(16) :: name='\psi_fra\psi_fra_'
        character(int_path+21) :: filename
        double precision :: x,psi_fra
        write(command,*) 'mkdir ',path,'\psi_fra'
        if(t==0) then
            call system(command)
        end if
        write(step,"(I5.5)") timestep/output_interval
        filename=path//name//step
        open(10,file=filename,status='replace')
        do i=1,num_ele
            x=dx*(i-1)
            psi_fra=Gc/2d0/L_0*((1d0-c(i))**2d0+L_0**2d0*((c(i+1)-c(i))/dx)**2d0)
            write(10,'(e24.12,1x,e24.12e4)') x,psi_fra
            x=dx*i
            psi_fra=Gc/2d0/L_0*((1d0-c(i+1))**2d0+L_0**2d0*((c(i+1)-c(i))/dx)**2d0)
            write(10,'(e24.12,1x,e24.12e4)') x,psi_fra
        end do
        close(10)
    end subroutine output_psi_fra
    !================================================================
    subroutine output_energy
        character(int_path+7) :: filename_1
        character(int_path+13) :: filename_2
        character(7) :: name_1='\energy'
        character(13) :: name_2='\energy_total'
        double precision :: kinetic,strain,fracture,total
        filename_1=path//name_1
        filename_2=path//name_2
        !変数初期化
        kinetic=0
        strain=0
        fracture=0
        !各種エネルギー密度を積分
        do i=1,num_ele
            kinetic=kinetic+density/2d0*(v(i+1)**2d0+v(i+1)*v(i)+v(i)**2d0)/3d0*dx
            strain=strain+(c(i+1)**2d0+c(i+1)*c(i)+c(i)**2d0)/3d0*(lamda/2d0+myu)*((u(i+1)-u(i))/dx)**2d0*dx
            fracture=fracture+Gc/2d0/L_0*((c(i+1)**2d0+c(i+1)*c(i)+c(i)**2d0-3d0*c(i+1)-3d0*c(i)+3d0)/3d0*dx&
            +L_0**2d0*((c(i+1)-c(i))/dx)**2d0*dx)
        end do
        total=kinetic+strain+fracture
        !出力
        if(t==0) then
            open(10,file=filename_1,status='replace')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') u(num_nod),kinetic,strain,fracture
            close(10)
            open(10,file=filename_2,status='replace')
            write(10,'(e24.12)') total
            close(10)
            !初期エネルギーの記録
            energy_initial=kinetic+strain+fracture
            energy_max=max(kinetic,strain,fracture)
        else
            open(10,file=filename_1,status='old',position='append')
            write(10,'(e24.12,1x,e24.12,1x,e24.12,1x,e24.12)') u(num_nod),kinetic,strain,fracture
            close(10)
            open(10,file=filename_2,status='old',position='append')
            write(10,'(e24.12)') total
            close(10)
            if(energy_max<max(kinetic,strain,fracture)) then
                energy_max=max(kinetic,strain,fracture)
            end if
        end if
    end subroutine output_energy
    !================================================================
    subroutine output_a_exa
        character(int_path+12) :: filename
        character(6) :: name='\a_exa'
        double precision :: t,x,a_exa
        filename=path//name
        open(10,file=filename,status='replace')
        do i=1,100001
            x=L_x/2d0
            t=4d0/10d0**9d0*(i-1)-L_x/2d0/v_d
            if(t>0) then
                a_exa=Amp*(Gauss(t,2)*sin(w*t)+2*Gauss(t,1)*w*cos(w*t)-Gauss(t,0)*w**2d0*sin(w*t))
                write(10,'(e24.12,1x,e24.12)') 4d0/10d0**9d0*(i-1),a_exa
            else
                write(10,'(e24.12,1x,e24.12)') 4d0/10d0**9d0*(i-1),0
            end if
        end do
        close(10)
    end subroutine output_a_exa
    !================================================================
    subroutine break_check
        if(count_3==0) then
            do i=1,num_nod
                if(c(i)<1d0/10d0**3d0) then
                    count_3=timestep/output_interval
                    write(*,*) 'Break!'
                    write(*,*) count_3
                    u_break_3=u(num_nod)
                    t_break_3=t
                end if
            end do
        end if
        if(count_4==0) then
            do i=1,num_nod
                if(c(i)<1d0/10d0**4d0) then
                    count_4=timestep/output_interval
                    write(*,*) 'Break!'
                    write(*,*) count_4
                    u_break_4=u(num_nod)
                    t_break_4=t
                end if
            end do
        end if
        if(count_5==0) then
            do i=1,num_nod
                if(c(i)<1d0/10d0**5d0) then
                    count_5=timestep/output_interval
                    write(*,*) 'Break!'
                    write(*,*) count_5
                    u_break_4=u(num_nod)
                    t_break_4=t
                end if
            end do
        end if
    end subroutine break_check
    !================================================================


    !================================================================
    !パラメーターをファイルに出力するサブルーチン
    !================================================================
    subroutine output_para
        character(11) :: name_p='\parameters'
        character(int_path+11) :: filename_p
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
        write(10,'(a,e24.12)') 'strain_c=',strain_c
        write(10,'(a,e24.12)') 'stress_c=',stress_c
        write(10,'(a,e24.12)') 'strain_weak=',(Gc_weak/3d0/L_0/(lamda+2d0*myu))**0.5d0
        write(10,'(a,e24.12)') 'stress_weak=',9d0/16d0*((lamda+2d0*myu)*Gc_weak/3d0/L_0)**0.5d0
        ! write(10,'(a,e24.12)') 'energy_max=',energy_max
        close(10)
    end subroutine
    !================================================================
    !ビデオを作るmatlabの.mファイルを作る
    !================================================================
    subroutine output_matlab
        character(79) :: filename_m
        character(5) :: step
        filename_m='C:\Users\tanaka\Documents\phase_field_1d_newest_code\phase_field_1d\makevideo.m'
        open(10,file=filename_m,status='replace')

        !\\\\\\\\\\\\\\\\\\\\\\\\\sig.mp4
        write(10,'(a)') "clear;"
        write(10,'(a)') "count=1;"
        write(10,'(a,i,a)') "num=",int(total_timestep/output_interval/10d0),";"
        write(10,'(a)') "fig=figure;"
        write(10,'(a)') "frames(num+1)=struct('cdata',[],'colormap',[]);"
        
        write(10,'(a)') "for i=0:num"

        write(10,'(a)') "if i>count*num/10"
        write(10,'(a)') "message=[sprintf('%2u0',count),'% is completed'];"
        write(10,'(a)') "disp(message);"
        write(10,'(a)') "count=count+1;"
        write(10,'(a)') "end"

        write(10,'(a)') "filenum=sprintf('%04u0',i);"
        write(10,'(a,a,a)') "filename=append('",path,"\sig\sig_',filenum);"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "[sigmax,index]=max(data(:,2));"
        write(10,'(a)') "plot(data(:,1),data(:,2),data(index,1),data(index,2),'pentagram');"
        write(10,'(a,e24.12,a)') "xlim([0 ",L_x,"]);"
        write(10,'(a,e24.12,1x,e24.12,a)') "ylim([",-1.1d0*stress_c,1.1d0*stress_c,"]);"
        write(10,'(a,e24.12,a)') "tim=sprintf('%.2e',",dt*output_interval*10d0,"*i);"
        write(10,'(a)') "time=append(tim,' [s]');"
        write(10,'(a)') "title(time);"
        write(10,'(a)') "drawnow;"
        write(10,'(a)') "frames(i+1)=getframe(fig);"
        write(10,'(a)') "end"
        write(10,'(a,a,a)') "video=VideoWriter('",path,"\sig.mp4','MPEG-4');"
        write(10,'(a)') "video.Quality=100;"
        write(10,'(a)') "video.FrameRate=160;"
        write(10,'(a)') "open(video);"
        write(10,'(a)') "writeVideo(video,frames);"
        write(10,'(a)') "close(video);"

        write(10,'(a)') "disp('sig.mp4 is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\c.mp4

        write(10,'(a)') "clear;"
        write(10,'(a)') "count=1;"
        write(10,'(a,i,a)') "num=",int(total_timestep/output_interval/10d0),";"
        write(10,'(a)') "fig=figure;"
        write(10,'(a)') "frames(num+1)=struct('cdata',[],'colormap',[]);"

        write(10,'(a)') "for i=0:num"
        
        write(10,'(a)') "if i>count*num/10"
        write(10,'(a)') "message=[sprintf('%2u0',count),'% is completed'];"
        write(10,'(a)') "disp(message);"
        write(10,'(a)') "count=count+1;"
        write(10,'(a)') "end"

        write(10,'(a)') "filenum=sprintf('%04u0',i);"
        write(10,'(a,a,a)') "filename=append('",path,"\c\c_',filenum);"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "xlim([0 ",L_x,"]);"
        write(10,'(a)') "ylim([-0.01 1.01]);"
        write(10,'(a,e24.12,a)') "tim=sprintf('%.2e',",dt*output_interval*10d0,"*i);"
        write(10,'(a)') "time=append(tim,' [s]');"
        write(10,'(a)') "title(time);"
        write(10,'(a)') "drawnow;"
        write(10,'(a)') "frames(i+1)=getframe(fig);"

        write(10,'(a)') "end"

        write(10,'(a,a,a)') "video=VideoWriter('",path,"\c.mp4','MPEG-4');"
        write(10,'(a)') "video.Quality=100;"
        write(10,'(a)') "video.FrameRate=160;"
        write(10,'(a)') "open(video);"
        write(10,'(a)') "writeVideo(video,frames);"
        write(10,'(a)') "close(video);"

        write(10,'(a)') "disp('c.mp4 is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\energy-u

        write(10,'(a)') "clear;"
        write(10,'(a)') "fig=figure;"

        write(10,'(a,a,a)') "filename=append('",path,"\energy');"
        write(10,'(a)') "data_1=load(filename);"
        write(10,'(a,a,a)') "filename=append('",path,"\energy_total');"
        write(10,'(a)') "data_2=load(filename);"
        write(10,'(a)') "plot(data_1(:,1),data_1(:,2),data_1(:,1),data_1(:,3),data_1(:,1),data_1(:,4),&
        data_1(:,1),data_2(:,1));"
        write(10,'(a,e24.12,a)') "xlim([0 ",u(num_nod),"]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "xline(",u_break_3,",'-.');"
        write(10,'(a)') "legend('kinetic','strain','fracture','total','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "print(fig,'-djpeg','",path,"\energy-u.jpg','-r600');"

        write(10,'(a)') "disp('energy-u.jpg is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\sig-u
        write(10,'(a)') "clear;"
        write(10,'(a)') "count=1;"
        write(10,'(a,i,a)') "n_3=",count_3,";"
        write(10,'(a,i,a)') "n_4=",count_4,";"
        write(10,'(a,i,a)') "n_5=",count_5,";"
        write(10,'(a,i,a)') "num=",int(total_timestep/output_interval),";"
        write(10,'(a,i,a)') "num_2=",1,";"
        write(10,'(a,i,a)') "num_3=",nint(num_ele*2d0/4d0),";"
        write(10,'(a,i,a)') "num_4=",num_ele,";"
        write(10,'(a,i,a)') "num_5=",nint(num_ele*2d0/4d0*3d0),";"
        write(10,'(a,i,a)') "num_6=",num_ele*2,";"
        write(10,'(a)') "sig=zeros(num+1,6);"
        write(10,'(a)') "c=zeros(num+1,6);"
        
        write(10,'(a)') "fig=figure;"
        
        write(10,'(a)') "for i=0:num"

        write(10,'(a)') "filenum=sprintf('%05u',i);"
        write(10,'(a,a,a)') "filename=append('",path,"\u\u_',filenum);"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a,i,a)') "sig(i+1,1)=data(",num_nod,",2);"

        write(10,'(a,a,a)') "filename=append('",path,"\sig\sig_',filenum);"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "sig(i+1,2)=data(num_2,2);"
        write(10,'(a)') "sig(i+1,3)=data(num_3,2);"
        write(10,'(a)') "sig(i+1,3)=(sig(i+1,3)+data(num_3+1,2))/2;"
        write(10,'(a)') "sig(i+1,4)=data(num_4,2);"
        write(10,'(a)') "sig(i+1,4)=(sig(i+1,4)+data(num_4+1,2))/2;"
        write(10,'(a)') "sig(i+1,5)=data(num_5,2);"
        write(10,'(a)') "sig(i+1,5)=(sig(i+1,5)+data(num_5+1,2))/2;"
        write(10,'(a)') "sig(i+1,6)=data(num_6,2);"

        write(10,'(a,a,a)') "filename=append('",path,"\c\c_',filenum);"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a,i,a)') "c(i+1,2)=data(",1,",2);"
        write(10,'(a,i,a)') "c(i+1,3)=data(",nint(num_nod/4d0),",2);"
        write(10,'(a,i,a)') "c(i+1,4)=data(",nint(num_nod/2d0),",2);"
        write(10,'(a,i,a)') "c(i+1,5)=data(",nint(num_nod*3d0/4d0),",2);"
        write(10,'(a,i,a)') "c(i+1,6)=data(",num_nod,",2);"

        write(10,'(a)') "end"

        write(10,'(a)') "plot(sig(1:n_3+1,1),sig(1:n_3+1,2),sig(1:n_3+1,1),sig(1:n_3+1,3)&
        ,sig(1:n_3+1,1),sig(1:n_3+1,4),sig(1:n_3+1,1),sig(1:n_3+1,5),sig(1:n_3+1,1),sig(1:n_3+1,6));"
        write(10,'(a)') "xlim([0 sig(n_3+1,1)]);"
        write(10,'(a,e24.12,a)') "ylim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(",stress_c,",'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_3
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_3+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_4+1,1),sig(1:n_4+1,2),sig(1:n_4+1,1),sig(1:n_4+1,3)&
        ,sig(1:n_4+1,1),sig(1:n_4+1,4),sig(1:n_4+1,1),sig(1:n_4+1,5),sig(1:n_4+1,1),sig(1:n_4+1,6));"
        write(10,'(a)') "xlim([0 sig(n_4+1,1)]);"
        write(10,'(a,e24.12,a)') "ylim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(",stress_c,",'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_4
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_4+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_5+1,1),sig(1:n_5+1,2),sig(1:n_5+1,1),sig(1:n_5+1,3)&
        ,sig(1:n_5+1,1),sig(1:n_5+1,4),sig(1:n_5+1,1),sig(1:n_5+1,5),sig(1:n_5+1,1),sig(1:n_5+1,6));"
        write(10,'(a)') "xlim([0 sig(n_5+1,1)]);"
        write(10,'(a,e24.12,a)') "ylim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(",stress_c,",'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_5
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_5+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\sig-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "disp('sig-u.jpg is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\c-u
        write(10,'(a)') "plot(sig(1:n_3+1,1),c(1:n_3+1,2),sig(1:n_3+1,1),c(1:n_3+1,3)&
        ,sig(1:n_3+1,1),c(1:n_3+1,4),sig(1:n_3+1,1),c(1:n_3+1,5),sig(1:n_3+1,1),c(1:n_3+1,6));"
        write(10,'(a)') "xlim([0 sig(n_3+1,1)]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_3
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_3+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_4+1,1),c(1:n_4+1,2),sig(1:n_4+1,1),c(1:n_4+1,3)&
        ,sig(1:n_4+1,1),c(1:n_4+1,4),sig(1:n_4+1,1),c(1:n_4+1,5),sig(1:n_4+1,1),c(1:n_4+1,6));"
        write(10,'(a)') "xlim([0 sig(n_4+1,1)]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_4
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_4+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_5+1,1),c(1:n_5+1,2),sig(1:n_5+1,1),c(1:n_5+1,3)&
        ,sig(1:n_5+1,1),c(1:n_5+1,4),sig(1:n_5+1,1),c(1:n_5+1,5),sig(1:n_5+1,1),c(1:n_5+1,6));"
        write(10,'(a)') "xlim([0 sig(n_5+1,1)]);"
        write(10,'(a,e24.12,a)') "xline(",strain_c*L_x,",'--');"
        write(10,'(a,e24.12,a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_5
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"
        write(10,'(a,e24.12,a)') "xlim([",strain_c*L_x*0.97," sig(n_5+1,1)]);"
        write(10,'(a,a,a)') "filename=append(name,'\c-u_",step,"_detail.jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "disp('c-u.jpg is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\c-sig
        write(10,'(a)') "plot(sig(1:n_3+1,2),c(1:n_3+1,2),sig(1:n_3+1,3),c(1:n_3+1,3)&
        ,sig(1:n_3+1,4),c(1:n_3+1,4),sig(1:n_3+1,5),c(1:n_3+1,5),sig(1:n_3+1,6),c(1:n_3+1,6));"
        write(10,'(a,e24.12,a)') "xlim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",stress_c,",'--');"
        write(10,'(a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_3
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-sig_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_4+1,2),c(1:n_4+1,2),sig(1:n_4+1,3),c(1:n_4+1,3)&
        ,sig(1:n_4+1,4),c(1:n_4+1,4),sig(1:n_4+1,5),c(1:n_4+1,5),sig(1:n_4+1,6),c(1:n_4+1,6));"
        write(10,'(a,e24.12,a)') "xlim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",stress_c,",'--');"
        write(10,'(a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_4
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-sig_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "plot(sig(1:n_5+1,2),c(1:n_5+1,2),sig(1:n_5+1,3),c(1:n_5+1,3)&
        ,sig(1:n_5+1,4),c(1:n_5+1,4),sig(1:n_5+1,5),c(1:n_5+1,5),sig(1:n_5+1,6),c(1:n_5+1,6));"
        write(10,'(a,e24.12,a)') "xlim([0 ",stress_c*1.1d0,"]);"
        write(10,'(a,e24.12,a)') "xline(",stress_c,",'--');"
        write(10,'(a)') "yline(0.75,'--');"
        write(10,'(a)') "legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(step,"(I5.5)") count_5
        write(10,'(a,a,a)') "name='",path,"';"
        write(10,'(a,a,a)') "filename=append(name,'\c-sig_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "disp('c-sig.jpg is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\c-x
        write(10,'(a,a,a)') "name='",path,"';"
        write(step,"(I5.5)") count_3
        write(10,'(a,a,a)') "filename=append(name,'\c\c_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\c_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(step,"(I5.5)") count_4
        write(10,'(a,a,a)') "filename=append(name,'\c\c_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\c_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(step,"(I5.5)") count_5
        write(10,'(a,a,a)') "filename=append(name,'\c\c_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\c_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "disp('c-x.jpg is created');"

        !\\\\\\\\\\\\\\\\\\\\\\\\\sig-x
        write(10,'(a,a,a)') "name='",path,"';"
        write(step,"(I5.5)") count_3
        write(10,'(a,a,a)') "filename=append(name,'\sig\sig_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\sig-x_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(step,"(I5.5)") count_4
        write(10,'(a,a,a)') "filename=append(name,'\sig\sig_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\sig-x_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(step,"(I5.5)") count_5
        write(10,'(a,a,a)') "filename=append(name,'\sig\sig_",step,"');"
        write(10,'(a)') "data=load(filename);"
        write(10,'(a)') "plot(data(:,1),data(:,2));"
        write(10,'(a,e24.12,a)') "L_x=",L_x,";"
        write(10,'(a)') "xlim([0 L_x]);"
        write(10,'(a)') "newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};"
        write(10,'(a)') "colororder(newcolors);"
        write(10,'(a)') "drawnow;"
        write(10,'(a,a,a)') "filename=append(name,'\sig-x_",step,".jpg');"
        write(10,'(a)') "print(fig,'-djpeg',filename,'-r600');"

        write(10,'(a)') "disp('sig-x.jpg is created');"

        close(10)
    end subroutine output_matlab
    !=================================================================
end program main