##utility that solves for the direction, speed, and time a projectile needs to hit a moving target
@tool class_name BallisticUtil extends RefCounted

const TIME_MAX: float = 20.0 ##the maximum time (s) the solver checks shots for- larger values find longer shots but are slower
const TIME_STEP: float = 0.01 ##the time step (s) between solves and the shot time for the first solve- larger values are faster but can skip over shots
const PRECISION_STEPS: int = 20 ##number of bisections to pin down exact crossing time, each one halves the error
const DAMP_MIN: float = 0.000001 ##drag below this is treated as no drag since the drag equations divide by it
const AXIS_MIN: float = 0.0001 ##the shortest a cross product can be and still be used as a rotation axis
const TIEBREAK: float = 0.001 ##cost added per second of flight time so shorter shots win when everything else is close

const WEIGHT_MISS:  float = 1.0  ##the miss distance (m) that costs 1.0- larger values make the miss matter less
const WEIGHT_TIME:  float = 1.0  ##the time error (s) from the time param that costs 1.0- larger values make it matter less
const WEIGHT_SPEED: float = 10.0 ##the speed error (m/s) from the speed param that costs 1.0- larger values make it matter less
const WEIGHT_DIR:   float = 0.5  ##the direction error (radians) from the dir param that costs 1.0- larger values make it matter less

enum Param {TIME, SPEED, DIR} ##the keys used for solve's params argument

#NOTE: wind is something that's set by Area3Ds so it shouldnt be a global constant- this is just
#here to show that it can be implemented
static var wind: Vector3 = Vector3.ZERO ##velocity (m/s) of the air, drag pulls the projectile towards it

#NOTE: this runs about TIME_MAX / TIME_STEP solves (plus bisections) per call, so don't call it every frame
##returns the best shot (direction, speed, and time) to hit a moving target, or null if there's no track or no valid shot
##pp and pv are the projectile's position and velocity (before launch kick)
##track stores the target's position, velocity, and acceleration, as of its last update
##damp is the linear RigidBody3D drag applied on the projectile
##grav is the gravitational acceleration to account for
##speed_range is the min (x) and max (y) speeds (m/s) the projectile can travel/be launched at
##params sets hints for what the solver tries to aim for- TIME (seconds), SPEED (m/s), and DIR (a Vector3 world direction,
##or a Vector2 offset in degrees from the direction towards the target, where x pitches up and y yaws right)
##thrust is a fixed motor acceleration vector on the projectile and steer is the motor acceleration (m/s^2) along the direction being solved for
static func solve(
		pp: Vector3,
		pv: Vector3,
		track: Track,
		damp: float,
		grav: Vector3 = ProjectSettings.get_setting("physics/3d/default_gravity") * ProjectSettings.get_setting("physics/3d/default_gravity_vector"),
		speed_range: Vector2 = Vector2(0.0, INF),
		params: Dictionary[Param, Variant] = {},
		thrust: Vector3 = Vector3.ZERO,
		steer: float = 0.0) -> Solution:
	if track == null: return null #without a target there's nothing to aim at
	#NOTE: all of this is explained in detail in the docs
	#wind goes thru drag (it pulls the projectile towards the wind's velocity) so it's an acceleration of wind * damp
	var a: Vector3 = grav + wind * damp + thrust
	var calc: Callable = func(t: float) -> Solution:
		#e replaces t and g replaces (1/2)t^2 because drag decays velocity and acceleration over time-
		#with no drag they equal the originals (the drag equations divide by damp so they can't be used)
		var e: float = t if damp < DAMP_MIN else (1.0 - exp(-damp * t)) / damp
		var g: float = t * t * 0.5 if damp < DAMP_MIN else (t - e) / damp
		#the gap the launch kick + steered motor has to bridge- where the target will be minus where
		#the projectile ends up if it was never launched
		var d: Vector3 = track.pos + track.vel * t + track.acc * (0.5 * t * t) - (pp + pv * e + a * g)
		var n: float = (d.length() - steer * g) / e #the speed needed to bridge the gap in time
		var s: float = clampf(n, speed_range.x, speed_range.y) #the speed that can actually be used
		#makes and returns the final solution- passing in the dir, speed, time, the speed (n) before
		#clamping, and how many meters the shot missed by (distance = velocity * time)
		return Solution.new(d.normalized(), s, t, n, absf(n - s) * e)
	
	#this array has all of the speeds where it matters to know the exact times where needed speed, n,
	#equals that value: the min, max, and speed param
	var spds: Array[float] = [speed_range.x]
	if speed_range.y < INF and speed_range.y != speed_range.x: spds.append(speed_range.y)
	if params.has(Param.SPEED): spds.append(params.get(Param.SPEED) as float)
	var sols: Array[Solution] #stores all of solutions from flight times and how far they miss
	#add the solution at the time parameter if one was provided, and clamp to the time limits
	if params.has(Param.TIME): sols.append(calc.call(clampf(params.get(Param.TIME) as float, TIME_STEP, TIME_MAX)))
	#add the first solution to the results
	var prev: Solution = calc.call(TIME_STEP)
	sols.append(prev)
	#since the 1st time step is already in the solutions array, start from the next one and step
	#all the way to the max (+1 to include it)
	for i: int in range(2, roundi(TIME_MAX / TIME_STEP) + 1):
		var cur: Solution = calc.call(i * TIME_STEP) #a solution is calculated for each time step and added to the results
		sols.append(cur)
		#for each speed level, check if the needed speed went from above the level to below it
		#or vice versa between the 2 samples. the 2 sides of the != differ only if 1 sample is above
		#the level and the other isnt- meaning the crossing time of needed and special speeds 
		#lies inside this time step window
		for v: float in spds:
			if (prev.need > v) != (cur.need > v):
				#we know the crossing is between the previous and current times, so we repeatedly
				#get the middle value and calculate the time there. we check the window again, and
				#if the crossing is in the upper window, we raise the low time to midpoint, and if
				#it's in the lower window, we drop the high time to the midpoint.
				var l: float = prev.time
				var h: float = cur.time
				for _i: int in PRECISION_STEPS:
					var m: float = (l + h) * 0.5
					if ((calc.call(m) as Solution).need > v) == (prev.need > v): l = m
					else: h = m
				sols.append(calc.call(h)) #after making precise, add the shot from the high window
		prev = cur #make the previous the current one for the next iteration
	
	#the direction param is turned into a world vector- if a vector3, assume it's already one, but
	#if a vector2, it's an offset (x pitches up and y yaws right, in degrees, with world up as +Y)
	#from the direction towards the target that gets turned into a world direction
	var dir: Vector3 = Vector3.ZERO
	if params.has(Param.DIR):
		if params.get(Param.DIR) is Vector3: dir = params.get(Param.DIR)
		else:
			var o: Vector2 = params.get(Param.DIR) as Vector2
			#using dir to target, it rotates only the yaw to the offset (- flips it for preference)
			var y: Vector3 = (track.pos - pp).normalized().rotated(Vector3.UP, deg_to_rad(-o.y))
			#axis to rotate around for the pitch tilt, this is horizontal to the yawed direction with +Y as up
			var r: Vector3 = y.cross(Vector3.UP)
			if r.length() < AXIS_MIN: r = Vector3.RIGHT #if the yawed direction is vertical, just use right as the tilting axis
			dir = y.rotated(r.normalized(), deg_to_rad(o.x)) #rotate on that tilting axis by the pitch offset
	var res: Solution = null #the best solution
	var rcst: float = INF #the weighted cost of the best solution
	for sol: Solution in sols: #for each solution,
		#cost begins at a small tiebreaker where shorter shots cost less, and the amount of the miss
		#is also weighed and added onto the cost
		var c: float = sol.time * TIEBREAK + pow(sol.miss / WEIGHT_MISS, 2.0)
		#we then add onto the cost by how far the shot is from the provided parameters if they
		#exist, and that error is multiplied by the weights. ^2 makes a big miss hurt more than
		#several small ones
		if params.has(Param.TIME):  c += pow((sol.time - (params.get(Param.TIME) as float)) / WEIGHT_TIME, 2.0)
		if params.has(Param.SPEED): c += pow((sol.speed - (params.get(Param.SPEED) as float)) / WEIGHT_SPEED, 2.0)
		if params.has(Param.DIR):   c += pow(sol.dir.angle_to(dir) / WEIGHT_DIR, 2.0)
		if c < rcst: #if the cost is less than the current best solution, update the current best solution and cost
			rcst = c
			res = sol
	return res #return the best solution

##holds the solution for a launcher or projectile to use
class Solution extends RefCounted:
	var dir: Vector3 ##the direction to point towards
	var speed: float ##the speed for the projectile to travel at
	var time: float ##the time it takes for the projectile to reach the target
	var need: float ##the raw needed speed, unclamped by speed_range
	var miss: float ##the distance in meters the projectile will miss by
	
	##stores the values of the solution
	func _init(_dir: Vector3, _speed: float, _time: float, _need: float, _miss: float) -> void:
		dir = _dir
		speed = _speed
		time = _time
		need = _need
		miss = _miss

##stores tracking data for a target
class Track extends RefCounted:
	var pos: Vector3 ##the last recorded position of the target
	var vel: Vector3 ##the velocity of the target
	var acc: Vector3 ##the acceleration of the target
	
	var _hist: Array[Vector3] ##the most recent positions recorded for the target
	var _cap: int ##the max length of the history
	
	##cap is the max amount of positions to keep in the history
	func _init(cap: int = 3) -> void: _cap = cap
	
	##updates position, velocity, and acceleration using the most recent position and the history
	##val is the new position and dt is the seconds since the last update (does nothing if it's 0 or less)
	func update(val: Vector3, dt: float) -> void:
		if dt <= 0.0: return #no time passed so velocity and acceleration can't be calculated
		#add to history, and if it's too long, remove the oldest position
		_hist.append(val)
		if _hist.size() > _cap: _hist.pop_front()
		var sz: int = _hist.size()
		pos = val #update latest position
		vel = (_hist[-1] - _hist[0]) / (dt * (sz - 1)) if sz > 1 else Vector3.ZERO #use first and last positions to get the velocity
		acc = ((_hist[-1] - _hist[-2]) - (_hist[1] - _hist[0])) / (dt * dt * (sz - 2)) if sz > 2 else Vector3.ZERO #use first and last velocities to get acceleration
