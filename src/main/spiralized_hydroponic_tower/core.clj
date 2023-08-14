(ns spiralized-hydroponic-tower.core
  (:refer-clojure :exclude [set])
  (:require
   [plexus.utils :as pu]
   [spiralized-hydroponic-tower.utils :as u]
   [spiralized-hydroponic-tower.math
    :refer [pi|2 pi|3 pi|4 pi|5 pi|6 two-pi
            cos acos sin asin tan atan pi
            TT T T|2 T|3 T|4 T|5 T|6]]
   [plexus.core :as paths
    :refer [forward hull left right up down backward defmodel translate spin
            set set-meta segment lookup-transform transform branch offset
            rotate frame save-transform add-ns extrude to points
            result difference union insert intersection loft export-models]]
   [clj-manifold3d.core :as m]))

(def cup-holder-inner-radius 53/2)
(def base-extra-radius 25)
(def base-height 205)
(def base-height-extra-tall 400)
(def tower-wall-thickness 1.0)
(def tower-segment-radius 57.8)
(def tower-segment-amplitude 12.7)

(def spacer-segment-height 208)
(def spacer-segment-offset (+ 2 1.2 base-extra-radius))
(def spacer-segment-angle (atan (/ spacer-segment-height spacer-segment-offset)))
(def spacer-segment-length (/ (abs spacer-segment-offset) (cos spacer-segment-angle)))

(def default-face-number 50)

(def fogger-base-radius (+ 47/2 1.6))
(def fogger-reservoir-height (+ 60 17))
(def fogger-reservoir-radius (+ fogger-base-radius 0.5))

(def top-joint--height 10)
(def top-joint-offset 0.5)
(def top-joint-angle (atan (/ 0.5 10)))

(def tube-outer-radius (+ 0.3 (/ (u/in->mm 5/16) 2)))
(def tube-inner-radius (- (/ (u/in->mm 5/16) 2) 1.6))

(defn housing-shape-pts
  ([circle-radius sin-wave-amplitude wall-thickness ratio]
   (housing-shape-pts circle-radius sin-wave-amplitude wall-thickness ratio 120))
  ([circle-radius sin-wave-amplitude wall-thickness ratio n-steps]
   (let [angle (* 2 pi)
         step-size (/ angle n-steps)]
     (for [step (range (inc n-steps))
           :let [x (* step step-size)
                 a (+ (+ circle-radius wall-thickness)
                      (* sin-wave-amplitude ratio (+ 1 (Math/sin (- (* 3 (- x (/ pi 2))) (* 3/3 pi)))))
                      (* sin-wave-amplitude (- 1 ratio) (+ 1 (Math/sin (* 3 (- x (/ pi 2)))))))]]
       [(* a (cos x))
        (* a (sin x))]))))

(defn housing-shape
  "Sum a circle function with four revolutions of a sin function."
  [circle-radius sin-wave-amplitude wall-thickness ratio]
  (m/cross-section (housing-shape-pts circle-radius sin-wave-amplitude wall-thickness ratio)))

(def tower-segment-contour
  (points
   :axes [:x :y]
   :meta-props [:segment]
   (frame :name :origin :meta {:segment :bottom})
   (rotate :x (- pi|2))
   (translate :x tower-segment-radius)
   (forward :length 160 :n-steps 55)
   (set-meta :segment :top)
   (u/curve-segment :offset 0.8
                    :height 3
                    :cs 4)
   (forward :length 5)
   (rotate :y top-joint-angle)
   (forward :length 10)))

(let [[bottom top] (partition-by (fn [x] (-> x meta :segment))
                                 tower-segment-contour)]
  (def tower-segment-bottom-contour bottom)
  (def tower-segment-top-contour top))

(def ^:export-model spiralized-tower-segment
  (extrude
   (result :name :spiralized-tower-segment
             :expr
             (union
              (difference (intersection :side-support :pod-segment-wall-outer) :net-pod-mask-composite)
              (difference :center-support :net-pod-mask-composite)
              (difference :pod-segment-wall-outer :pod-segment-wall-inner :net-pod-mask-composite)
              (intersection :pod-segment-wall-outer
                            (difference :net-pod-wall :net-pod-mask-composite))))

   (result :name :net-pod-wall
             :expr
             (union
              (for [i (range 3)]
                (rotate :z (* 2/3 i pi)
                        (union (hull (union :net-pod-wall-body)
                                     (translate :x 80 (union :net-pod-wall-body)))
                               (hull (union :net-pod-wall-2 :net-pod-wall-body)
                                     (translate :x 20 (union :net-pod-wall-2))))))))

   (result :name :net-pod-mask-composite
           :expr
           (union
            (for [i (range 3)]
              (rotate {:z (* 2/3 i pi)}
                      (union :net-pod-mask
                             (hull :net-pod-mask-2 (translate {:x 30} (union :net-pod-mask-2))))))))

   (frame :name :origin)

   (branch
    :from :origin
    :with []
    (frame :cross-section (m/circle cup-holder-inner-radius 70)
           :name :net-pod-mask)
    (frame :cross-section (m/circle (+ cup-holder-inner-radius 1.0 6) 30)
           :name :net-pod-wall-body)
    (translate :x (- tower-segment-radius 9.15))
    (forward :length 2)
    (offset :delta 6 :to [:net-pod-mask])
    (forward :length 3)
    (branch
       :from :net-pod-wall-body
       :with []
       (frame :cross-section (m/circle (+ cup-holder-inner-radius 1.0 6) 70)
              :name :net-pod-wall-2)
       (frame :cross-section (m/circle (+ 6 cup-holder-inner-radius) 70)
              :name :net-pod-mask-2)
       (forward :length 20)
       (right :angle pi|4 :curve-radius 150)
       (forward :length 45)))

   (branch
    :from :origin
    :with []
    (frame :name :center-support
           :cross-section (m/square 3 24 true))
    (rotate :z pi|3)
    (for [i (range 3)]
      (branch
       :from :origin :with [:center-support]
       (rotate :z (* i 2/3 pi))
       (translate :x 20)
       (forward :length 2))))

   (branch
    :from :origin
    :with []
    (frame :name :side-support
           :cross-section (m/square 3 30 true))
    (for [i (range 3)]
      (branch
       :from :origin
       (rotate :z (* i 2/3 pi))
       (translate :x (- 50))
       (segment
        (for [sign [- +]]
          (branch
           :from :origin
           (translate :y (sign 5))
           (rotate :z (sign (- pi|4)))
           (forward :length 4)))))))

   (branch
    :from :origin
    :with []
    (frame :name :pod-segment-wall-outer)
    (frame :name :pod-segment-wall-inner)
    (->> (concat
          (for [[radius amplitude] [[(- tower-segment-radius 1) 12.5]]]
            [(set :cross-section (housing-shape radius amplitude 0 0) :to [:pod-segment-wall-outer])
             (set :cross-section (housing-shape radius amplitude (- tower-wall-thickness) 0) :to [:pod-segment-wall-inner])
             (forward :length 0.05)])
          (for [[i [radius z]] (map-indexed list tower-segment-bottom-contour) ]
            (let [ratio (/ i (dec (count tower-segment-bottom-contour)))]
              [(set :cross-section (housing-shape radius tower-segment-amplitude 0 ratio) :to [:pod-segment-wall-outer])
               (set :cross-section (housing-shape radius tower-segment-amplitude (- tower-wall-thickness) ratio) :to [:pod-segment-wall-inner])
               (translate :z (+ 5 z) :global? true)
               (forward :length 0.05)]))
          (for [[_ [radius z]] (map-indexed list tower-segment-top-contour)]
            (let [ratio 1]
              [(set :cross-section (housing-shape radius tower-segment-amplitude 0 ratio) :to [:pod-segment-wall-outer])
               (set :cross-section (housing-shape radius tower-segment-amplitude (- tower-wall-thickness) ratio) :to [:pod-segment-wall-inner])
               (translate :z (+ 5 z) :global? true)
               (forward :length 0.05)])))
         (partition 2 1)
         (map (fn [[l r]]
                (hull l r))))
    (to
     :models [:pod-segment-wall-inner]
     (hull
      (translate :z -0.25)
      (forward :length 0.3))))))

(def ^:export-model fogger-reservoir
  (extrude
   (result :name :fogger-reservoir
           :expr (difference :fogger-reservoir-body :fogger-reservoir-mask))
   (frame :name :origin :fn 100)

   (branch
    :from :origin
    :with []
    (frame :name :fogger-reservoir-body
           :cross-section (m/circle fogger-reservoir-radius 100))
    (forward :length 1.8)
    (frame :name :fogger-reservoir-mask
           :cross-section (m/circle (- fogger-reservoir-radius 1.6) 100))
    (forward :length (- fogger-reservoir-height 12 35))
    (u/right-forward-left-segment
     :bottom-radius fogger-reservoir-radius
     :height 14
     :offset 7
     :curve-height 1.5
     :to [:fogger-reservoir-body])
    (u/right-forward-left-segment
     :bottom-radius (- fogger-reservoir-radius 1.6)
     :height 14
     :offset 7
     :curve-offset -1.6
     :curve-height 1.5
     :to [:fogger-reservoir-mask])
    (forward :length 3)
    (forward :length 0.1 :to [:fogger-reservoir-mask])
    (for [_ (range 3)]
      (hull
       (set :cross-section (m/circle (+ 7.8 fogger-reservoir-radius) 4) :to [:fogger-reservoir-body])
       (set :cross-section (m/circle (+ 6.9 fogger-reservoir-radius) 4) :to [:fogger-reservoir-mask])
       (forward :length 1)
       (translate :z 12)
       (set :cross-section (m/circle (+ 7 fogger-reservoir-radius) 120) :to [:fogger-reservoir-body])
       (set :cross-section (m/circle (+ 6.2 fogger-reservoir-radius) 120) :to [:fogger-reservoir-mask])
       (forward :length 1 :fn default-face-number)))
    (forward :length 5))))

(def ^:export-model spiralized-fogger-segment
  (extrude
   (result :name :fogger-reservoir
           :expr (union (difference :base :fog-holes)
                        (difference :pod-segment-wall-outer :pod-segment-wall-inner)))
   (frame :name :origin)

   (branch
    :from :origin
    (:forms spiralized-tower-segment))

   (branch
    :from :origin
    :with []
    (frame :name :base
           :cross-section (housing-shape (- tower-segment-radius 1) 12.7 0 0))
    (hull
     :to [:base]
     (forward :length 0.01)
     (forward :length 1.19)))

   (branch
    :from :origin
    (frame :cross-section (m/circle 6 20) :name :fog-holes)
    (rotate :z pi|3)
    (for [r [(- tower-segment-radius 10) (- tower-segment-radius 25)]
          [x y] (housing-shape-pts r 12.5 0 1 6)]
      (branch
       :from :origin
       :with [:fog-holes]
       (translate :x x :y y)
       (forward :length 20))))))

(def spacer-segment-contour
  (points
   :axes [:x :y]
   :meta-props [:segment]
   (frame :name :origin :meta {:segment :bottom})
   (rotate :x (- pi|2))
   (translate :x (+ tower-segment-radius base-extra-radius 2))
   (rotate :y (- (- pi|2 spacer-segment-angle)))
   (forward :length spacer-segment-length :n-steps 25)
   (forward :length 0.1)
   (set-meta :segment :top)
   (rotate :y (- pi|2 spacer-segment-angle))
   (u/curve-segment :offset 2 :height 6)
   (rotate :y top-joint-angle)
   (forward :length 10)))

(defn make-tower-base [base-height]
  (let [base-segment-contour (points
                              :axes [:x :y]
                              :meta-props [:segment]
                              (frame :name :origin :meta {:segment :bottom})
                              (rotate :x (- pi|2))
                              (translate :x (+ tower-segment-radius base-extra-radius 2))
                              (forward :length 10)
                              (u/curve-segment :offset (- 1.35)
                                               :height base-height)
                              (set-meta :segment :top)
                              (u/curve-segment :offset 2 :height 6 :cs 5)
                              (rotate :y top-joint-angle)
                              (forward :length 10))
        [base-tower-bottom-contour
         base-tower-top-contour] (partition-by (fn [x] (-> x meta :segment))
                                               base-segment-contour)]
    (extrude
     (result :name :spiralized-tower-base
             :expr (difference :base-wall-outer :base-wall-inner))
     (frame :name :base-wall-outer :cross-section (housing-shape (+ base-extra-radius tower-segment-radius 2) tower-segment-amplitude 0 0))
     (hull
      (forward :length 1.1)
      (forward :length 0.1))
     (frame :name :base-wall-inner)
     (->> (concat
           (for [[i [radius z]] (map-indexed list base-tower-bottom-contour)]
             (let [ratio (/ i (dec (count base-tower-bottom-contour)))]
               [(set :cross-section (housing-shape radius tower-segment-amplitude (* ratio 0.8) ratio) :to [:base-wall-outer])
                (set :cross-section (housing-shape radius tower-segment-amplitude (- (* ratio 0.8) tower-wall-thickness) ratio) :to [:base-wall-inner])
                (translate :z (+ 1.2 z) :global? true)
                (forward :length 0.01)]))
           (for [[_ [radius z]] (map-indexed list base-tower-top-contour)]
             (let [ratio 1]
               [(set :cross-section (housing-shape radius tower-segment-amplitude (* ratio 0.8) ratio) :to [:base-wall-outer])
                (set :cross-section (housing-shape radius tower-segment-amplitude (- (* ratio 0.8) tower-wall-thickness) ratio) :to [:base-wall-inner])
                (translate :z (+ 1.2 z) :global? true)
                (forward :length 0.01)])))
          (partition 2 1)
          (map (fn [[b t]]
                 (hull b t)))))))

(def ^:export-model spiralized-tower-base
  (make-tower-base base-height))

(def ^:export-model spiralized-tower-base-extra-tall
  (make-tower-base base-height-extra-tall))

(def ^:export-model spiralized-tower-lid
  (extrude
   (result :name :spiralized-tower-lid
           :expr (difference (union :lid-wall-outer :tube-body)
                             :lid-net-pod-mask :tube-mask))
   (frame :name :origin :fn 40)

   (branch
    :from :origin
    :with []
    (frame :name :tube-mask
           :cross-section (m/union
                           (for [i (range 3)]
                             (-> (m/hull (m/circle tube-inner-radius 100)
                                         (-> (m/circle tube-inner-radius 100)
                                             (m/translate [30 0])))
                                 (m/rotate (* i 2/3 T))))))
    (translate :z 1.2)
    (forward :length 3)
    (frame :name :tube-body
           :cross-section (m/circle tube-outer-radius 100))
    (set :cross-section (m/circle tube-inner-radius 100) :to [:tube-mask])
    (forward :length 4)
    (for [_ (range 3)]
      (hull
       (forward :length 1)
       (translate :z 3)
       (offset :delta -0.3 :to [:tube-body])
       (forward :length 0.1)
       (offset :delta 0.3 :to [:tube-body]))))

   (branch
    :from :origin
    :with []
    (frame :name :lid-wall-outer)
    (frame :name :lid-wall-inner)
    (->> (concat
          (for [amplitude [12.9 12.5]]
            [(set :cross-section (housing-shape tower-segment-radius amplitude 0 0) :to [:lid-wall-outer])
             (set :cross-section (housing-shape tower-segment-radius amplitude (- tower-wall-thickness) 0) :to [:lid-wall-inner])
             (forward :length (if (= amplitude 12.5) 8 0.1))])
          #_(for [[i [radius z]] (map-indexed list tower-segment-bottom-contour)]
              (let [ratio (/ i (dec (count tower-segment-bottom-contour)))]
                [(set :cross-section (housing-shape radius tower-segment-amplitude 0 ratio) :to [:lid-wall-outer])
                 (set :cross-section (housing-shape radius tower-segment-amplitude (- tower-wall-thickness) ratio) :to [:lid-wall-inner])
                 (translate :z (+ 7 z) :global? true)
                 (forward :length 0.01)])))
         (partition 2 1)
         (map (fn [[b t]]
                (hull b t)))))

   (branch
    :from :origin
    :with []
    (frame :cross-section (m/square 200 1.6 true)
           :name :base-supports)
    (for [i (range 3)]
      (branch
       :from :base-supports
       (rotate :z (* i 2/3 pi))
       (forward :length 4))))

   (branch
    :from :origin
    :with []
    (frame :cross-section (m/circle cup-holder-inner-radius)
           :name :lid-net-pod-mask)
    (for [i (range 3)]
      (branch
       :from :origin
       (rotate :z (* i 2/3 pi))
       (translate :x (- tower-segment-radius 9))
       (forward :length 10))))

   (branch
    :from :origin
    :with []
    (frame :name :tower-lid-top :cross-section (housing-shape tower-segment-radius tower-segment-amplitude 0 0))
    (hull
     (forward :length 1.2)
     (forward :length 0.1)))))

(def spiralized-tower-fogger-segment
  (extrude
   (result :name :spiralized-tower-fogger-segment
           :expr (union
                  :reservoir-connection-walls
                  :reservoir-bottom-connections
                  (difference :outer-wall
                              :inner-wall)
                  (difference :reservoir-outer-wall
                              :reservoir-inner-wall)))

   (result :name :reservoir-bottom-connections
             :expr (intersection :outer-wall
                                 (union
                                  (for [i (range 3)]
                                    (rotate :z (* i 2/3 pi) (union :reservoir-bottom-connection))))))


   (result :name :reservoir-connection-walls
           :expr (union
                  (for [i (range 3)]
                    (rotate :z (* i 2/3 pi) (union :reservoir-connection-wall)))))

   (frame :name :origin :fn 40)

   (let [w (- (+ 51 2.2) fogger-reservoir-radius)]
     (branch
      :from :origin
      :with []
      (frame :name :reservoir-connection-wall
             :cross-section (m/square w 1.6 true))
      (rotate :z pi|3)
      (translate :x (+ (/ w 2) fogger-reservoir-radius))
      (hull
       (forward :length 3)
       (translate :x (- (/ w 2)) :z (- fogger-reservoir-height 2))
       (set :cross-section (m/square 1.6 1.6 true))
       (forward :length 0.1))))

   (branch
    :from :origin
    :with []
    (frame :name :reservoir-bottom-connection
           :cross-section (m/square 50 15 true))
    (rotate :z pi|3)
    (translate :x 30)
    (forward :length 1.2))

   (branch
    :from :origin
    :with []
    (frame :name :outer-wall :cross-section (housing-shape tower-segment-radius tower-segment-amplitude (- tower-wall-thickness) 0))
    (frame :name :inner-wall :cross-section (housing-shape tower-segment-radius tower-segment-amplitude (- (* 2 tower-wall-thickness)) 0))
    (hull
     (forward :length 2.9)
     (forward :length 0.1)))

   (branch
    :from :origin
    :with []
    (frame :name :reservoir-outer-wall :cross-section (m/circle fogger-reservoir-radius))
    (forward :length 1.2)
    (frame :name :reservoir-inner-wall :cross-section (m/circle (- fogger-reservoir-radius 1.0)))
    (forward :length fogger-reservoir-height))))

(def mini-fan
  (extrude
   (result :name :mini-fan :expr (difference :fan-body :fan-bolt-holes :fan-blade-mask))
   (frame :cross-section (m/square 25 25 true) :name :fan-body :fn 70)
   (frame :cross-section (m/union
                          (for [x [10 -10]
                                y [10 -10]]
                            (-> (m/circle 3/2)
                                (m/translate [x y]))))
         :name :fan-bolt-holes)
   (frame :cross-section (m/circle (dec 25/2)) :name :fan-blade-mask)
   (forward :length 10)))

(def spacer-segment-contour
  (points
   :axes [:x :y]
   :meta-props [:segment]
   (frame :name :origin :meta {:segment :bottom})
   (rotate :x (- pi|2))
   (translate :x (+ tower-segment-radius base-extra-radius 2))
   (rotate :y (- (- pi|2 spacer-segment-angle)))
   (forward :length spacer-segment-length :n-steps 25)
   (rotate :y (- pi|2 spacer-segment-angle))
   (set-meta :segment :top)
   (u/curve-segment :offset 2
                    :height 6
                    :cs 4)
   (rotate :y top-joint-angle)
   (forward :length 10)))

(let [[bottom top] (partition-by (fn [x] (-> x meta :segment))
                                 spacer-segment-contour)]
  (def spacer-bottom-contour bottom)
  (def spacer-top-contour top))

(def ^:export-model spiralized-tower-spacer-segment
  (extrude
   (result :name :spiralized-tower-spacer-segment
           :expr (union
                  #_(difference (intersection :cable-slot-body :spacer-wall-outer) :cable-slot-mask)
                  #_(difference (intersection :base-supports :spacer-wall-outer) :cable-slot-mask :nozzle-mask)
                  (difference (intersection :spacer-wall-outer :spacer-net-pod-wall) :spacer-net-pod-mask-composite)
                  :base-latice
                  #_(difference (intersection :spacer-wall-outer :base-supports) :spacer-net-pod-mask-composite)
                  (difference :spacer-wall-outer :spacer-wall-inner :spacer-net-pod-mask-composite)))

   (result :name :base-latice
           :expr (difference :base-latice-body :spacer-net-pod-wall :base-holes))

   (result :name :spacer-net-pod-wall
           :expr
           (hull (union :spacer-net-pod-wall-body :spacer-net-pod-wall-2)
                 (translate :x 30 (union :spacer-net-pod-wall-body :spacer-net-pod-wall-2))))

   (result :name :spacer-net-pod-mask-composite
           :expr
           (union :spacer-net-pod-mask
                  (hull (union :spacer-net-pod-mask-2)
                        (translate :x 50 (union :spacer-net-pod-mask-2))) ))

   (frame :name :origin :fn 80)

   (branch
    :from :origin
    :with []
    (frame :name :base-supports
           :cross-section (m/union
                           (m/circle 22)
                           (m/union
                            (for [i (range 3)]
                              (-> (m/square (+ 50 tower-segment-radius) 20 true)
                                  (m/translate [(/ (+ 10 tower-segment-radius) 2) 0])
                                  (m/rotate (+ T|3 (* i 2/3 T))))))) )
    (frame :name :nozzle-mask
           :cross-section (m/circle 19))
    (frame :name :tube-mask
           :cross-section (m/circle (/ (u/in->mm 1/2) 2)))
    (translate :x -10 :to [:tube-mask])
    (forward :length 1.2)
    (forward :length 10 :to [:tube-mask]))

   (branch
    :from :origin
    :with []
    (frame :cross-section (m/circle (+ 18 cup-holder-inner-radius))
           :name :spacer-net-pod-mask)
    (frame :cross-section (m/circle (+ 18 cup-holder-inner-radius 1.0 6))
           :name :spacer-net-pod-wall-body)
    (for [i (range 1)]
      (branch
       :from :origin
       :with [:spacer-net-pod-mask :spacer-net-pod-wall-body]
       (rotate :z (* i 2/3 pi))
       (translate :x (- (+ 7 tower-segment-radius) 7))
       (forward :length 2)
       (offset :delta 6 :to [:spacer-net-pod-mask])
       (forward :length 3)
       (branch
        :from :spacer-net-pod-wall-body
        :with []
        (frame :cross-section (m/circle (+ 18 cup-holder-inner-radius 1.0 6))
               :name :spacer-net-pod-wall-2)
        (frame :cross-section (m/circle (+ 18 6 cup-holder-inner-radius))
               :name :spacer-net-pod-mask-2)
        (forward :length 10)
        (right :angle pi|6 :curve-radius 200)
        (forward :length 60)))))

   (branch
    :from :origin
    :with []
    (frame :name :cable-slot-body :cross-section (m/square 9 12 true))
    (frame :name :cable-slot-mask :cross-section (m/square 6 9 true))
    (translate :x (+ tower-segment-radius 1 (* 2 tower-segment-amplitude)))
    (forward :length 12)
    (right :angle (/ pi 6) :curve-radius 9/2)
    (forward :length 10))

   (branch
    :from :origin
    :with []
    (frame :name :base-latice-body
           :cross-section (housing-shape (- (+ tower-segment-radius base-extra-radius 1.8)
                                            0.8)
                                         (- tower-segment-amplitude 0.0) 0 0))
    (hull
     (forward :length 0.01)
     (translate :z 1.58)
     (forward :length 0.01)))

   (branch
    :from :origin
    :with []
    (frame :name :base-holes :cross-section (m/circle (/ (u/in->mm 1/2) 2)))
    (segment
     (for [[x y] [[-10 0]
                  [-10 70]
                  [-50 70]
                  [-50 0]
                  [-50 -70]
                  [-10 -70]]]
       (branch
        :from :base-holes
        (translate :x x :y y)
        (forward :length 4)))))

   (branch
    :from :origin
    :with []
    (frame :name :spacer-wall-outer)
    (frame :name :spacer-wall-inner)
    (->> (concat
          (for [[amplitude length] [[(- tower-segment-amplitude 0.4) 5] [(+ 0.3 tower-segment-amplitude) 0.1]]]
              [(set :cross-section (housing-shape (+ tower-segment-radius base-extra-radius 2) amplitude 0 0) :to [:spacer-wall-outer])
               (set :cross-section (housing-shape (+ tower-segment-radius base-extra-radius 2) amplitude (- tower-wall-thickness) 0) :to [:spacer-wall-inner])
               (forward :length length)])
          (for [[i [radius z]] (map-indexed list spacer-bottom-contour)]
            (let [ratio (/ i (dec (count spacer-bottom-contour)))]
              [(set :cross-section (housing-shape radius tower-segment-amplitude 0 ratio) :to [:spacer-wall-outer])
               (set :cross-section (housing-shape radius tower-segment-amplitude (- tower-wall-thickness) ratio) :to [:spacer-wall-inner])
               (translate :z (+ 5 z) :global? true)
               (forward :length 0.01)]))
          (for [[_ [radius z]] (map-indexed list spacer-top-contour)]
            (let [ratio 1]
              [(set :cross-section (housing-shape radius tower-segment-amplitude 0 ratio) :to [:spacer-wall-outer])
               (set :cross-section (housing-shape radius tower-segment-amplitude (- tower-wall-thickness) ratio) :to [:spacer-wall-inner])
               (translate :z (+ 5 z) :global? true)
               (forward :length 0.01)])))
         (partition 2 1)
         (map (fn [[b t]]
                  (hull b t)))))))

(def ^:export-model big-fan-mount
  (extrude
   (result :name :big-fan-mount :expr (difference :mount-body :mount-mask :cable-slot-mask))

   (frame :name :origin :fn 50)

   (branch
    :from :origin
    :with []
    (frame :name :cable-slot-mask
           :cross-section #_(m/circle 11)
           (m/union
              (for [i (range 2)]
                (-> (m/square 20 0.5 true)
                    (m/rotate (* i 1/2 T))))))
    (translate :x 10)
    (forward :length 1.2)
    (set :cross-section (m/circle 13) :to [:cable-slot-mask])
    (forward :length 20))

   (branch
    :from :origin
    :with []
    (frame :cross-section (m/circle (+ 21 cup-holder-inner-radius)) :name :mount-body :fn 60)
    (frame :cross-section (m/square 25 25 true) :name :mount-mask)
    (translate :x (- 20) :to [:mount-mask])
    (forward :length 2)
    (offset :delta -3 :to [:mount-body])
    (forward :length 3)
    (set :cross-section (m/square (+ 25 1.8) (+ 25 1.8) true)
         :to [:mount-body])
    (translate :x -20 :to [:mount-body])
    (forward :length 4)
    (hull
     (forward :length 0.1)
     (translate :z 1.8)
     (offset :delta -1.0 :to [:mount-mask])
     (forward :length 0.1)))))

(def fan-mount
  (extrude
   (result :name :fan-mount
           :expr (difference :fan-mount-body :fan-mount-mask :fan-bolt-holes :fan-blade-mask))
   (frame :name :fan-mount-body
          :cross-section (m/square (+ 25 (* 2 1.2)) (+ 25 (* 2 1.2)) true))
   (branch :from :fan-mount-body :with [] (:forms mini-fan))
   (forward :length 1)
   (frame :name :fan-mount-mask
          :cross-section (m/square 25 25 true))
   (forward :length 15)
   (hull
    (forward :length 1)
    (translate :z 15)
    (set :cross-section (m/circle (+ 19 0.5) 100) :to [:fan-mount-body])
    (set :cross-section (m/circle (- (+ 19 0.5) 0.5) 100) :to [:fan-mount-mask])
    (forward :length 0.1))
   (for [[model curve-offset radius]
         [[:fan-mount-body 0 (+ 19 0.5)]
          [:fan-mount-mask -0.5 19]]]
     [(u/curved-cylinder-offset
       :bottom-radius radius
       :offset -0.5
       :height 1.5
       :curve-offset curve-offset
       :cs 100
       :to [model])
      (forward :length 0.8 :to [model])
      (u/curved-cylinder-offset
       :bottom-radius (- radius 0.5)
       :offset 0.5
       :height 1.5
       :curve-offset (- curve-offset)
       :cs 100
       :to [model])])
   (forward :length 1/2 :to [:fan-mount-mask])))

(def ^:export-model full-tower
  (extrude
   (result :name :full-tower
           :expr (union :spiralized-tower-base
                        (->> (union :spiralized-tower-spacer-segment)
                             (translate :z 230)
                             (rotate :z pi|3))
                        (for [i (range 5)]
                          (->> (union :spiralized-tower-segment)
                               (translate :z (+ 230 220 (* i 160)))
                               (rotate :z (* i pi|3))))
                        (->> :fogger-reservoir
                             (translate {:z (+ 230 220 (* 5 160))})
                             (rotate {:z (* 5 pi|3)}))
                        (->> :spiralized-tower-lid
                             (translate {:z (+ 230 220 (* 6 160) 10)})
                             (rotate {:z (* 6 pi|3)}))))
   (insert :extrusion spiralized-tower-base)
   (insert :extrusion spiralized-tower-spacer-segment)
   (insert :extrusion spiralized-tower-segment)
   (insert :extrusion spiralized-fogger-segment)
   (insert :extrusion spiralized-tower-lid)))

(comment

  (time (doseq [model (export-models *ns* "glb")] (deref model)))
  (time (doseq [model (export-models *ns* "stl")] (deref model)))

  )
