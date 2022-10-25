program findgravity

    implicit none

    real(8), dimension(6) :: top_heights = (/0.527, 0.527, 0.526, 0.526, 0.526, 0.526/)
    real(8), dimension(6) :: bottom_heights = (/0.291, 0.416, 0.364, 0.182, 0.488, 0.383/)
    real(8), dimension(6) :: top_time = (/0.03658, 0.03277, 0.03865, 0.03592, 0.0388, 0.03728/)
    real(8), dimension(6) :: bottom_time = (/0.1901, 0.12381, 0.15579, 0.23347, 0.07741, 0.14381/)

    integer :: n

    do n = 1, 6
        call airdrag(top_heights(n), bottom_heights(n), top_time(n), bottom_time(n))
    end do

end program findgravity

subroutine airdrag(top_height, bottom_height, top_time, bottom_time)
    ! Uses gradient descent to find gravitation acceleration that minimizes error.
    implicit none

    real(8) :: loss, get_acceleration

    real(8) :: top_height ! Distance to the Top of the Apparatus
    real(8) :: bottom_height ! Distance to the Bottom of the Apparatus
    real(8) :: top_time ! Time to Unblock the Top Laser
    real(8) :: bottom_time ! Time to Unblock the Bottom Laser

    real(8), parameter :: delta_g = 0.1
    real(8) :: learning_rate = 0.1
    real(8), parameter :: minimum_ratio = 1.0 + 1e-7

    real(8) :: gradient, gravitational_acceleration, current_loss, old_loss
    real(8) :: negative_dir, positive_dir, pos_top, pos_bottom, neg_top, neg_bottom

    integer(8) :: n

    gravitational_acceleration = -1.0 ! Arbitrary Initial Condition, such that the Simulation isn't Unstable.
    old_loss = 2.0
    current_loss = 1.0
    n = 0

    do while (old_loss / current_loss > minimum_ratio)
        old_loss = current_loss
        n = n + 1
        call simulate_fall(gravitational_acceleration - delta_g, top_height, bottom_height, neg_top, neg_bottom)
        call simulate_fall(gravitational_acceleration + delta_g, top_height, bottom_height, pos_top, pos_bottom)
        negative_dir = loss(top_time, neg_top, bottom_time, neg_bottom)
        positive_dir = loss(top_time, pos_top, bottom_time, pos_bottom)
        current_loss = 0.5 * (negative_dir + positive_dir)
        write(100, *) n, current_loss
        gradient = (positive_dir - negative_dir) / (2.0 * delta_g)
        gravitational_acceleration = gravitational_acceleration - learning_rate * gradient

    end do

    write(*,*) gravitational_acceleration

    return

end subroutine airdrag

real(8) function loss(exp_top, act_top, exp_bottom, act_bottom)

    implicit none

    real(8) :: exp_top, act_top, exp_bottom, act_bottom ! Variables Passed
    real(8) :: top_diff, bottom_diff ! Internal Variables for Faster Computation Time

    top_diff = exp_top - act_top
    bottom_diff = exp_bottom - act_bottom
    loss = top_diff * top_diff + bottom_diff * bottom_diff

    return

end function loss

subroutine simulate_fall(gravitational_acceleration, top_height, bottom_height, output_top, output_bottom)

    implicit none

    real(8) :: gravitational_acceleration, output_top, output_bottom
    real(8) :: top_height, bottom_height

    real(8), parameter :: pi = 4.0*atan(1.0)
    real(8), parameter :: height_max = 100000000. ! Max Height to Stop Program (in case of wrong way gravity)

    real(8), parameter :: delta_time = 0.001 ! Time Increment in Seconds
    real(8), parameter :: ball_mass = 0.01 ! Ball Mass in Kilograms
    real(8), parameter :: ball_diameter = 0.025 ! Ball Diameter in Meters
    real(8), parameter :: drag_coefficient = 0.47 ! Drag Coefficient of Sphere
    real(8), parameter :: air_density = 1.225 ! Density of Air in kg/m^3
    real(8), parameter :: distance_start = 0.009 ! Distance to Laser in Meters
    real(8), parameter :: starting_velocity = 0.0 ! Starting Velocity in Meters per Second

    real(8), parameter :: reference_area = pi * ball_diameter * ball_diameter * 0.25
    real(8) :: start
    real(8) :: activate_time
    real(8) :: middle
    real(8) :: finish

    real(8) :: last_pos_y, pos_y, last_vel_y, vel_y, acc_y, force_y, time, pass_top, pass_bottom
    real(8) :: a, b ! variables to do some linear interpolation for extra accuracy for basically free
    logical :: checkpointed, time_started

    start = top_height + distance_start + ball_diameter
    activate_time = top_height + 0.9374449248709414 * ball_diameter ! Can be lower
    middle = top_height + 0.06255507512905856 * ball_diameter ! Can be higher
    finish = bottom_height + 0.06255507512905856 * ball_diameter ! Can be higher

    pos_y = start
    vel_y = starting_velocity
    acc_y = gravitational_acceleration
    time = 0.0

    checkpointed = .FALSE.
    time_started = .FALSE.

    do while ((pos_y > finish) .AND. (pos_y < height_max))

        last_pos_y = pos_y
        last_vel_y = vel_y
        time = time + delta_time
        call rk4(pos_y, vel_y, delta_time, air_density, drag_coefficient, reference_area, ball_mass,&
                gravitational_acceleration)
        if (pos_y < activate_time .AND. (.NOT. time_started)) then
            a = (activate_time - pos_y) / (last_pos_y - pos_y)
            b = 1.0-a
            time = 0.0
            time_started = .TRUE.
            pos_y = activate_time
            vel_y = a * last_vel_y + b * vel_y
        end if
        if (pos_y < middle .AND. (.NOT. checkpointed)) then
            a = (middle - pos_y) / (last_pos_y - pos_y)
            b = 1.0-a
            pass_top = a * (time - delta_time) + b * time
            checkpointed = .TRUE.
        end if

    end do

    a = (finish - pos_y) / (last_pos_y - pos_y)
    b = 1.0-a
    pass_bottom = a * (time - delta_time) + b * time

    output_top = pass_top
    output_bottom = pass_bottom

    return

end subroutine simulate_fall

subroutine rk4(pos, vel, time_step, d, dc, ra, m, g)

    implicit none

    real(8) :: get_acceleration
    real(8), parameter :: sixth = 1.0/6.0

    real(8) :: d,  dc, ra, m, g ! Values to pass to get_acceleration
    real(8) :: pos, vel, time_step ! Values used for RK4
    real(8) :: pos_temp, vel_temp, acc_temp ! Variables used internally
    real(8), dimension(2) :: k1, k2, k3, k4

    pos_temp = pos
    vel_temp = vel
    acc_temp = get_acceleration(d, vel_temp, dc, ra, m, g)

    k1(1) = vel_temp
    k1(2) = acc_temp

    pos_temp = pos + 0.5 * time_step * k1(1)
    vel_temp = vel + 0.5 * time_step * k1(2)
    acc_temp = get_acceleration(d, vel_temp, dc, ra, m, g)

    k2(1) = vel_temp
    k2(2) = acc_temp

    pos_temp = pos + 0.5 * time_step * k2(1)
    vel_temp = vel + 0.5 * time_step * k2(2)
    acc_temp = get_acceleration(d, vel_temp, dc, ra, m, g)

    k3(1) = vel_temp
    k3(2) = acc_temp

    pos_temp = pos + time_step * k3(1)
    vel_temp = vel + time_step * k3(2)
    acc_temp = get_acceleration(d, vel_temp, dc, ra, m, g)

    k4(1) = vel_temp
    k4(2) = acc_temp

    pos = pos + sixth * time_step * (k1(1) + 2.0 * k2(1) + 2.0 * k3(1) + k4(1))
    vel = vel + sixth * time_step * (k1(2) + 2.0 * k2(2) + 2.0 * k3(2) + k4(2))

    return

end subroutine rk4

real(8) function get_acceleration(density, velocity, drag_coefficient, reference_area, mass, gravity)

    implicit none

    real(8) :: density, velocity, drag_coefficient, reference_area, mass, gravity, acceleration
    acceleration = (0.5 * density * velocity * velocity * drag_coefficient * reference_area) / mass
    get_acceleration = acceleration + gravity
    return

end function get_acceleration
