(ns kalman-sim.core
  (:use [incanter.core]
        [incanter.stats]
        [incanter.charts]
        [kalman-sim.simulator]))

(def size 100)

(defn init-state []
  (let [x (- (rand (* 2 size)) size)
        y (- (rand (* 2 size)) size)
        theta (- (rand (* 2 Math/PI)) Math/PI)]
    (dosync
     (ref-set state [x y theta]))))

(def Q (matrix [[1e-3 0 0] 
                [0 1e-3 0]
                [0 0 1e-4]])) ; Cov in executing commands
(def R (matrix [[0.2 0]
                [0 0.05]])) ; Cov in measurment

(defn wrap-angle [theta]
  (if (> theta Math/PI)
    (- theta (* 2 Math/PI))
    (if (<= theta (- Math/PI) )
      (+ theta (* 2 Math/PI))
      theta)))

(defn g
  "The function modeling state transitions\n
  u: command vector\n
  x: state vector"
  [u x]
  (let [[x y theta] x
        [v omega] u
        x (+ x  (* v (Math/cos theta) delta-t))
        y (+ y  (* v (Math/sin theta) delta-t))
        theta (wrap-angle (+ theta (* omega delta-t)))]
    [x y theta]))

(defn G
  "The Jacobian of g\n
  u: command vector\n
  x: state vector"
  [u x]
  (let [[x y theta] x
        [v omega] u]
    (matrix [[1 0 (* -1 v (Math/sin theta) delta-t)]
             [0 1 (* v (Math/cos theta) delta-t)]
             [0 0 1]])))

(defn h
  "The function modeling measurment\n
  x: state vector"
  [x]
  (let [[x y theta] x
        d (Math/sqrt (+ (* x x) (* y y)))
        phi (- (* Math/PI 1.5) theta (Math/atan2 y x))]
    [d phi]))

(defn H
  "The Jacobian of h\n
  x: state vector"
  [state]  (let [[x y theta] state
                 l2 (let [s (+ (* x x) (* y y))]
                      (if (== 0 s) 1e-9 s))]
             (matrix [[(/ x (Math/sqrt l2)) (/ y (Math/sqrt l2)) 0]
                      [(/ y l2)             (- (/ x l2))         -1]])))

(defn planner [mu sigma]
  [0 0])

(defn dist [[x y theta]]
  (Math/sqrt (+ (* x x) (* y y))))

(defn sample-data [log-entry]
  (sample-mvn 100 :mean (:mu log-entry) :sigma (:sigma log-entry)))

(defn plot-belief []
  (let [entries (take-nth 100 @history)
        plot (box-plot (map dist (sample-data (first entries)))
                       :x-label "Number of Measurments"
                       :y-label "Distance from Landmark"
                       :title "Belief")] 
    (doseq [entry (rest entries)]
      (add-box-plot plot (map dist (sample-data entry))))
    (view plot)))

(defn plot-true-state []
  (let [entries (take-nth 100 @history)
        plot (xy-plot (range (count entries)) (map #(dist (:state %1)) entries)
                   :x-label "Number of Measurments"
                   :y-label "Distance from Landmark"
                   :title "True State")]
    (set-x-label plot "Number of Measurments")
    (set-y-label plot "Distance from Landmark")
    (set-title plot "True State" )
    (view plot)))

(defn -main
  "Runs the simulation, outputs a graph when the trace of the
  covariance matrix is small enough."
  [& args]
  (run-simulation! g h G H Q R [(/ size 2) (/ size 2) 0] (matrix [[size 0 0]
                                                                 [0 size 0]
                                                                 [0 0 Math/PI]])
                   planner (init-state))
  (plot-belief)
  (plot-true-state))