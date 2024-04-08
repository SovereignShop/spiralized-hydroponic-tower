(ns spiralized-hydroponic-tower.utils
  (:require
   [spiralized-hydroponic-tower.math
    :refer [vec-add vec-dot vec-sub vec-magnitude
            vec-normalize cross-product-z vec-scale
            angle-between-vectors vec-length
            cos acos sin asin tan atan sqrt
            pi pi|2 pi|3 pi|4 pi|5 pi|6
            TT T T|2 T|3]] 
   [clj-manifold3d.core :as m]
   [plexus.transforms :as tf] 
   [sicmutils.calculus.derivative :as cd]
   [sicmutils.env :as e]
   [plexus.core :as p]
   [plexus.triangles :as triangles]))

(def epsilon 1e-6)

(defn near-zero? [x]
  (< (Math/abs x) epsilon))

(defn collinear?
  [p1 p2 p3]
  (near-zero? (cross-product-z (vec-sub p2 p1) (vec-sub p3 p1))))

(defmacro pts [& forms]
  `(p/points
    :axes [:x :z]
    (p/frame :name :origin)
    ~@forms))

(defn in->mm [inches]
  (* inches 25.4))

(defn ft->mm [ft]
  (* ft 12 25.4))

(defn mm->in [mm]
  (/ mm 25.4))

(defn frame-2d->point [frame]
  (-> frame (.getColumn 2) vec))

(defn get-neighbors [pts idx]
  (let [n (count pts)
        p0 (if (zero? idx) (dec n) (dec idx))
        p2 (if (= idx (dec n)) 0 (inc idx))]
    [(nth pts p0) (nth pts idx) (nth pts p2)]))

(defn circle-pts [r angle res]
  (let [angle (min (* 2 pi) angle)
        angle-inc (/ angle res)
        frame (m/frame-2d)
        pts (for [i (range (if (>= angle (- (* 2 pi) epsilon)) res (inc res)))]
              (-> frame
                  (m/rotate (* i angle-inc))
                  (m/translate [r 0])
                  (.getColumn 2)
                  (vec)))]
    (vec pts)))

(defn generate-curve-points
  [start-frame radius curve-angle n-segments]
  (let [angle-inc (/ curve-angle (dec n-segments))
        op (if (pos? curve-angle) 1 -1)
        pivot-frame (m/translate start-frame [0 (* op (- radius))])]
    (map (fn [i]
           (-> pivot-frame
               (m/rotate (* i angle-inc))
               (m/translate [0 (* op radius)])))
         (range n-segments))))

(defn curve-mid-point
  [start-frame radius curve-angle n-segments]
  (let [rot-angle (/ curve-angle 2)
        op (if (pos? curve-angle) 1 -1)
        pivot-frame (m/translate start-frame [0 (* op (- radius))])]
    (-> pivot-frame
        (m/rotate rot-angle)
        (m/translate [0 (* op radius)])
        (.getColumn 2)
        (vec))))

(defn abc->ABC
  "side,side,side -> angle,angle,angle using Law of Cosines."
  [a b c]
  (let [C (Math/acos (/ (- (Math/pow c 2)
                           (Math/pow a 2)
                           (Math/pow b 2))
                        (- (* 2 a b))))
        B (Math/acos (/ (- (Math/pow b 2)
                           (Math/pow a 2)
                           (Math/pow c 2))
                        (- (* 2 a c))))
        A (- Math/PI (+ C B))]
    [A B C]))

(defn print-defs [m]
  (doseq [[k v] m]
    (println (format "(def %s %s)" (str (name k)) (str v)))))

(defn three-point-circle [p1 p2 p3]
  (let [v1 (vec-sub p2 p1)
        v2 (vec-sub p3 p1)

        a (vec-length v1)
        b (vec-length (vec-sub p3 p2))
        c (vec-length v2)
        [A B C] (abc->ABC a b c)
        radius (/ a (* 2 (sin A)))
        D (acos (/ a (* 2 radius)))


        direction (cross-product-z v1 v2)

        center
        (-> (m/frame-2d (vec-normalize (vec-sub p2 p1)) p1)
            (m/rotate (if (pos? direction) D  (- D)))
            (m/translate [radius 0])
            (frame-2d->point))]
    [radius center direction]))

(defn three-point-arc [p1 p2 p3 n-segments]
  (let [[radius center direction] (three-point-circle p1 p2 p3)
        v1 (vec-sub p1 center)
        v2 (vec-sub p3 center)
        direction2 (cross-product-z v1 v2)
        curve-angle (cond->> (angle-between-vectors
                              (vec-sub p1 center)
                              (vec-sub p3 center))
                      false (- (* 2 pi)))
        angle-increment (/ curve-angle n-segments)
        start-dir (vec-normalize (vec-sub p1 center))
        frame (m/frame-2d start-dir center)]
    (for [i (range (inc n-segments))]
      (-> frame
          (m/rotate (* i (if (neg? direction) -1 1) angle-increment))
          (m/translate [radius 0])))))

(comment

  (let [[p1 p2 p3 :as pts] [[10 0] [0 10] [-10 0]] #_[[79.89636993408203 -17.5]
                                                    [14.585479736328125 1.646345140215999E-7]
                                                    [79.89636993408203 17.5]]
      frames
      (three-point-arc
       p1 p2 p3 20)]
  (m/difference
   #_(-> (m/circle radius 40)
       (m/translate center))
   (m/union (for [pt (map frame-2d->point frames)]
              (-> (m/circle 2) (m/translate pt))))))

(let [[p1 p2 p3 :as pts] [[79.89636993408203 -17.5]
                          [14.585479736328125 1.646345140215999E-7]
                          [79.89636993408203 17.5]]
      frames
      (three-point-arc
       p1 p2 p3 20)]
  (m/difference
   #_(-> (m/circle radius  40)
         (m/translate center))
   (m/union (for [pt (map frame-2d->point frames)]
              (-> (m/circle 3) (m/translate pt))))
   (m/union (for [pt pts]
            (-> (m/circle 3) (m/translate pt))))))





  (m/cross-section
   (map frame-2d->point
        (three-point-arc [10 0] [0 10] [-10 0] 30)))

  (let [[p1 p2 p3]
        (for [rot [(double 0) pi|4 (+ pi|4 pi|2)]]
          (-> (m/frame-2d 1)
              (m/translate [0 20])
              (m/rotate rot)
              (m/translate [10 0])
              (.getColumn 2)
              (vec)))]
    (m/union
     (m/cross-section
      (map frame-2d->point
           (three-point-arc p1 p2 p3 20)))
     (m/union
      (for [p [p1 p2 p3]]
        (-> (m/circle 1/2 30)
            (m/translate p))))))


  )

(defn round-edge
  [p1 p2 p3 radius n-segments]
  (if (collinear? p1 p2 p3)
    (take n-segments (repeat p2))
    (let [vdist (vec-sub p2 p1)
          v1 (vec-normalize vdist)
          v2 (vec-normalize (vec-sub p3 p2))
          direction (cross-product-z v1 v2)
          curve-angle (cond-> (angle-between-vectors v1 v2)
                        (not (pos? direction)) -)
          start-pos (vec-add p1 (vec-scale v1 (- (vec-length vdist)
                                                 (* radius (tan (/ curve-angle 2))))))
          start-frame (m/frame-2d v1 start-pos)
          ;; mid-point (curve-mid-point start-frame radius curve-angle n-segments)
          curve-points (generate-curve-points
                        start-frame radius
                        curve-angle
                        n-segments)]
      curve-points)))

(defn round-polygon
  [points curve-radius n-segments]
  (for [[p1 p2 p3] (partition 3 1 (conj (vec points) (nth points 0) (nth points 1)))
        pt (round-edge p1 p2 p3 curve-radius n-segments)]
    pt))

(defn calc-right-forward-left-segment
  ([delta height curve-height curve-offset]
   (calc-right-forward-left-segment delta height curve-height curve-offset nil))
  ([delta height curve-height curve-offset n-forward-steps]
   (let [mid [(/ delta 2) (/ height 2)]
         hv [0 curve-height]
         ->mid (vec-sub mid hv)
         A (acos (/ (vec-dot hv ->mid)
                    (* (vec-magnitude hv) (vec-magnitude ->mid))))
         B (/ (- pi A) 2)
         C (- pi|2 B)

         radius (/ curve-height (tan C))
         angle (* 2 C)
         length (* 2 (- (vec-magnitude ->mid) (vec-magnitude hv)))]
     [(if (neg? delta)
        (p/left :angle angle :curve-radius (- radius curve-offset))
        (p/right :angle angle :curve-radius (- radius curve-offset)))
      (p/forward :length length :n-steps n-forward-steps)
      (if (neg? delta)
        (p/right :angle angle :curve-radius (+ radius curve-offset))
        (p/left :angle angle :curve-radius (+ radius curve-offset)))])))

(defn right-forward-left-segment
  [& {:keys [bottom-radius delta height curve-height curve-offset cs to]
      :or {curve-offset 0
           cs           100}}]
  (p/loft
   :to to
   (for [[r dz] (->> (p/points
                      :axes [:x :y]
                      (p/frame :name :origin)
                      (p/rotate :x (- pi|2))
                      (p/translate :x bottom-radius)
                      (calc-right-forward-left-segment delta height curve-height curve-offset))
                     (partition 2 1)
                     (map (fn [[[_ z0] [r z1]]]
                            [r (- z1 z0)])))]
     [(p/set :cross-section (m/circle r cs) :to to)
      (p/forward :length dz :to to)])))

(defn curve-segment
  [& {:keys [delta
             height
             curve-offset
             cs]
      :or {curve-offset 0
           cs 10}}]
  (if (zero? delta)
    (p/forward :length height)
    (let [h (/ height 2)
          w (/ delta 2)
          b (atan (/ h w))
          c (atan (/ w h))
          d (- b c)
          r (/ h (cos d))
          a (cond-> (asin (/ h r))
              (pos? delta) -)]
      [(p/left :angle a :curve-radius (+ r curve-offset) :cs cs)
       (p/right :angle a :curve-radius (- r curve-offset) :cs cs)])))

(defn curve-segment-2 [& {:keys [bottom-radius delta curve-offset height to cs vs angle]
                          :or {curve-offset 0
                               cs 100
                               vs 30
                               angle 360}}]
  (p/loft
   :to to
   (for [[r dz] (->> (p/points
                      :axes [:x :y]
                      (p/frame :name :origin)
                      (p/rotate :x (- pi|2))
                      (p/translate :x bottom-radius)
                      (curve-segment :delta delta :height height :curve-offset curve-offset :cs vs))
                     (partition 2 1)
                     (map (fn [[[_ z0] [r z1]]]
                            [r (- z1 z0)])))]
     [(p/set :cross-section (if (< angle 360)
                              (m/cross-section
                               (cons [0 0] (circle-pts r (/ angle 57.29578) cs)))
                              (m/circle r cs))
             :to to)
      (p/forward :length dz :to to)])))


(defn cos-wave-points
  ([width height]
   (cos-wave-points width height 1))
  ([width height n-waves]
   (cos-wave-points width height n-waves 0))
  ([width height n-waves delta]
   (let [steps 20
         step-size (+ 0.02 (/ (* 2 pi) steps))
         length (* n-waves 2 pi)
         dcos (cd/D e/cos)
         f (fn [x]
             (let [dv [(dcos x) -1]
                   ;; normal perp vector to tangent line, scaled by delta
                   [dx dy] (if (or (pos? delta) (neg? delta))
                             (e/* (- delta) (e// dv (e/magnitude dv)))
                             [0 0])]
               [(+ (- dx) (* (/ width 2) (/ x pi)))
                (+ (/ height 2) (- dy) (* (/ height 2) (Math/cos x)))]))]
     (conj (vec
            (for [x (range (- pi) (- length pi) step-size)]
              (f x)))
           (f (- length pi))))))

(defn make-gear
  [& {:keys [n-teeth
             tooth-height
             tooth-width]}]
  (let [c (* n-teeth tooth-width)
        r (/ c (* 2 pi))
        m tf/identity-tf
        pts (map (fn [[x y]]
                   [(+ x (/ tooth-width 2)) y])
                 (cos-wave-points tooth-width tooth-height n-teeth))
        n-pts (count pts)
        a (/ (* 2.00 pi) n-pts)]
    (m/cross-section
     (for [[i [_ y]] (map list (range n-pts) (butlast pts))]
       (-> m
           (tf/rotate :z (* i a))
           (tf/go-forward (- r y) :x)
           (tf/translation-vector)
           (subvec 0 2))))))

(defn ovol
  ([rx ry]
   (ovol rx ry 30))
  ([rx ry n-steps]
   (m/cross-section
    (for [x (range n-steps)]
      (let [d (* x (/ (* 2 Math/PI) n-steps))]
        [(* rx (Math/cos d))
         (* ry (Math/sin d))])))))

(defn show-pts [pts r]
  (m/union
   (for [pos pts]
     (-> (m/circle r 10)
         (m/translate pos)
         (m/extrude 1/2)))))

(defn arc-segment [side-length curve-radius props]
  (let [angle (triangles/abc->A side-length curve-radius curve-radius)
        r (- (/ Math/PI 2) (/ (- Math/PI angle) 2))]
    [(p/rotate :y (- r))
     (p/right :angle angle :curve-radius curve-radius :props props :cs 2)
     (p/rotate :y (- r))]))

(defn thread
  ([section pitch radius steps-per-revolution n-revolutions]
   (m/loft (m/union (for [i (range n-revolutions)]
                      (-> section
                          (m/translate [0 (- (* i pitch))]))))
           (for [step (range (inc steps-per-revolution))]
             (-> (m/frame 1)
                 (m/rotate [0 0 (* (/ (* 2 Math/PI) steps-per-revolution) step)])
                 (m/translate [radius 0 (* step (/ pitch steps-per-revolution))])
                 (m/rotate [(- (/ Math/PI 2)) 0 0])))))
  ([section pitch radius steps-per-revolution n-revolutions taper-slope]
   (m/loft section
           (for [rev (range n-revolutions)
                 step (range (inc steps-per-revolution))]
             (let [z (+ (* rev pitch) (* step (/ pitch steps-per-revolution)))
                   taper-offset (- z (* z taper-slope))]
               (-> (m/frame 1)
                   (m/rotate [0 0 (* (/ (* 2 Math/PI) steps-per-revolution) step)])
                   (m/translate [(- radius taper-offset) 0 z])
                   (m/rotate [(- (/ Math/PI 2)) 0 0])))))))

(defn screw-by-pitch
  [& {:keys [d0 length dr pitch flat taper]
      :or {flat 1}}]
  (let [r (/ d0 2)]
    (cond-> (m/difference
             (-> (m/union
                  (m/cylinder (+ length pitch) (- r dr) (* (or taper 1) (- r dr)) 100)
                  (if taper
                    (thread (ovol dr (/ pitch 2) 4) pitch (- r dr) 100 (+ 1 (/ length pitch)) taper)
                    (thread (ovol dr (/ pitch 2) 4) pitch (- r dr) 100 (+ 1 (/ length pitch)))))
                 (m/translate [0 0 (- (/ pitch 2))]))
             (-> (m/cylinder (* 2 pitch) (+ r 5))
                 (m/translate [0 0 length]))
             (-> (m/cylinder (* 2 pitch) (+ r 5))
                 (m/translate [0 0 (- (* 2 pitch))])))
      (< flat 1) (m/difference (m/union (let [m (-> (m/square (+ d0 5) d0 true)
                                                    (m/translate [0 (+ (/ d0 2) (/ d0 2) (- (/ (* (- 1 flat) d0) 2)))])
                                                    (m/extrude (+ length 5)))]
                                          (m/union m (m/mirror m [0 1 0]))))))))


(defn slice-model
  [model start end]
  (m/difference
   (m/translate model [0 0 (- start)])
   (m/union
    (-> (m/cube 1000 1000 1000 true)
        (m/translate [0 0 (+ 500 (- end start))]))
    (-> (m/cube 1000 1000 1000 true)
        (m/translate [0 0 (- 500)])))))

(defn triangle-pts
  [a b c]
  (let [[A B C] (triangles/abc->ABC a b c)]
    (next
     (p/points
      :axes [:x :z]
      (p/frame :name :origin)
      (p/rotate :y pi|2)
      (p/forward :length a)
      (p/rotate :y (- (- pi C)))
      (p/forward :length b)
      (p/rotate :y (- (- pi A)))
      (p/forward :length c)))))


(defn rounded-triangle-pts [a b c cr]
  (let [[A B C] (triangles/abc->ABC a b c)
        f (fn corner-offset [angle]
             (let [ang (/ angle 2)
                   length (/ cr (tan ang))]
              length))]
    (next
     (p/points
      :axes [:x :z]
      (p/frame :name :origin)
      (p/rotate :y pi|2)
      (p/forward :length (- a (f C)))
      (p/left :angle (- pi C) :curve-radius cr)
      (p/forward :length (- b (f C) (f A)))
      (p/left :angle (- pi A) :curve-radius cr)
      (p/forward :length (- c (f A) (f B)))
      (p/left :angle (- pi B) :curve-radius cr)))))

(defn triangle
  [a b c]
  (m/cross-section (triangle-pts a b c)))

(defn radians->degrees [radians]
  (* radians 57.29578))
