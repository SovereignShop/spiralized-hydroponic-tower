(ns spiralized-hydroponic-tower.utils
  (:require
   [spiralized-hydroponic-tower.math
    :refer [vec-add vec-dot vec-sub vec-magnitude
            cos acos sin asin tan atan
            pi pi|2 pi|3 pi|4 pi|5 pi|6
            TT T T|2 T|3]]
   [sicmutils.env :as e]
   [sicmutils.calculus.derivative :as cd]
   [sicmutils.generic :refer [dot-product magnitude negate cross-product]]
   [clj-manifold3d.core :as m]
   [plexus.transforms :as tf]
   [plexus.core :as p]
   [plexus.triangles :as triangles]))

(defn in->mm [inches]
  (* inches 25.4))

(defn ft->mm [ft]
  (* ft 12 25.4))

(defn mm->in [mm]
  (/ mm 25.4))

(defn circle-pts [r angle res]
  (let [pts (for [i (range 0 res)]
              (-> tf/identity-tf
                  (tf/rotate :z (* i (/ angle (dec res))))
                  (tf/go-forward r :x)
                  (tf/translation-vector)
                  (subvec 0 2)))]
    pts))

(defn calc-right-forward-left-segment
  ([offset height curve-height curve-offset]
   (calc-right-forward-left-segment offset height curve-height curve-offset nil))
  ([offset height curve-height curve-offset n-forward-steps]
   (let [mid [(/ offset 2) (/ height 2)]
         hv [0 curve-height]
         ->mid (vec-sub mid hv)
         A (acos (/ (vec-dot hv ->mid)
                    (* (vec-magnitude hv) (vec-magnitude ->mid))))
         B (/ (- pi A) 2)
         C (- pi|2 B)

         radius (/ curve-height (tan C))
         angle (* 2 C)
         length (* 2 (- (vec-magnitude ->mid) (vec-magnitude hv)))]
     [((if (neg? offset) p/left p/right) :angle angle :curve-radius (- radius curve-offset))
      (p/forward :length length :n-steps n-forward-steps)
      ((if (neg? offset) p/right p/left) :angle angle :curve-radius (+ radius curve-offset))])))

(defn right-forward-left-segment
  [& {:keys [bottom-radius offset height curve-height curve-offset cs to]
      :or {curve-offset 0
           cs           100}}]
  (p/loft
   :to to
   (for [[r dz] (->> (p/points
                      :axes [:x :y]
                      (p/frame :name :origin)
                      (p/rotate :x (- pi|2))
                      (p/translate :x bottom-radius)
                      (calc-right-forward-left-segment offset height curve-height curve-offset))
                     (partition 2 1)
                     (map (fn [[[_ z0] [r z1]]]
                            [r (- z1 z0)])))]
     [(p/set :cross-section (m/circle r cs) :to to)
      (p/forward :length dz :to to)])))

(defn curve-segment
  [& {:keys [offset
             height
             curve-offset
             cs]
      :or {curve-offset 0
           cs 10}}]
  (if (zero? offset)
    (p/forward :length height)
    (let [h (/ height 2)
          w (/ offset 2)
          b (atan (/ h w))
          c (atan (/ w h))
          d (- b c)
          r (/ h (cos d))
          a (cond-> (asin (/ h r))
              (pos? offset) -)]
      [(p/left :angle a :curve-radius (+ r curve-offset) :cs cs)
       (p/right :angle a :curve-radius (- r curve-offset) :cs cs)])))

(defn curved-cylinder-offset [& {:keys [bottom-radius offset curve-offset height to cs vs angle]
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
                      (curve-segment :offset offset :height height :curve-offset curve-offset :cs vs))
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
  ([width height n-waves offset]
   (let [steps 20
         step-size (+ 0.02 (/ (* 2 pi) steps))
         length (* n-waves 2 pi)
         dcos (cd/D e/cos)
         f (fn [x]
             (let [dv [(dcos x) -1]
                   ;; normal perp vector to tangent line, scaled by offset
                   [dx dy] (if (or (pos? offset) (neg? offset))
                             (e/* (- offset) (e// dv (e/magnitude dv)))
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
