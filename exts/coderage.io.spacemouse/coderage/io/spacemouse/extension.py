from datetime import datetime, timedelta
import math
import platform

import carb
import omni.ext
import omni.kit.viewport.window
import omni.kit.app
import omni.kit.viewport
import omni.ui.scene
import omni.ui as ui
import omni.usd
# import omni.kit.viewport.utility as vp_util
from omni.kit.viewport.utility import get_active_viewport
from omni.kit.viewport.utility.camera_state import ViewportCameraState as VpCamera
from pxr import Usd, UsdGeom, Gf, Sdf
if platform.system() == 'Windows':
    omni.kit.pipapi.install("pywinusb")
# import pywinusb
import spacenavigator

UPDATE_TIME_MILLIS = 10
DEFAULT_STEP = 50          # Ideally we'd get this from Camera Speed @TODO
DEFAULT_ROTATION = 1
DEBUG = True

# Spacemouse supports six degrees of freedom. By default, these are mapped to the ViewPort camera as so:
#  * x: tracking (trucking) left and right
#  * y: dollying forwards and backwards (move the *camera's* Omniverse-z axis)
#  * z: pedestal lift/lower
#  * roll: rolling the camera body on its z axis (the line-of-sight axis)
#  * pitch: tilting the camera up and down (rotating camera body on its horizontal axis)
#  * yaw: panning left/right (rotating camera body on its vertical axis)

class CoderageIoSpacemouseExtension(omni.ext.IExt):

    def on_startup(self, ext_id):
        self._count = 0
        self.__previous_time = None
        self.__previous_state = None
        viewport = get_active_viewport()
        self.__camera_path = str(viewport.camera_path)
        self._MovementValue = 1.0
        self._RotationValue = 1.0
        self._MovementScalar = DEFAULT_STEP
        self._RotationScalar = DEFAULT_ROTATION
        print("[coderage.io.spacemouse] coderage io spacemouse startup")

        # Angles
        def RollLeft_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": self._RotationValue, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def PanRight_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": self._RotationValue, "buttons": [0,0]}
            )
            self.update_state(state)

        def TiltDown_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": self._RotationValue, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def RollRight_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": -self._RotationValue, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def PanLeft_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": -self._RotationValue, "buttons": [0,0]}
            )
            self.update_state(state)

        def TiltUp_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": -self._RotationValue, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        # Movements
        def Up_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": self._MovementValue, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def Forward_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": self._MovementValue, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def Down_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": -self._MovementValue, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def Left_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": -self._MovementValue, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def Back_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": -self._MovementValue, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def Right_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": self._MovementValue, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state)

        def XAxisUp_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": self._MovementValue, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        def YAxisUp_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": self._MovementValue, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        def ZAxisUp_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": self._MovementValue, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        def XAxisDown_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": -self._MovementValue, "y": 0.0, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        def YAxisDown_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": -self._MovementValue, "z": 0.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        def ZAxisDown_Click():
            state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
                **{"t": 255.0, "x": 0.0, "y": 0.0, "z": -self._MovementValue, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
            )
            self.update_state(state, True)

        self._window = ui.Window("Spacemouse debug", width=300, height=300)
        with self._window.frame:
            with ui.VStack():

                # add an IntSlider for translate Strength
                ui.Label("Camera Rotation Amount")
                self._rotationSlider = ui.IntSlider(min = 1, max = 90, step=5)                
                self._rotationSlider.model.set_value(self._RotationScalar)
                # self._rotationValue = 5
                self._rotationSlider.model.add_value_changed_fn(self._onrotation_value_changed)

                with ui.HStack():                    
                    rollLeftButton = ui.Button("Roll L", clicked_fn=RollLeft_Click)
                    panLeftButton = ui.Button("Pan L", clicked_fn=PanLeft_Click)
                    tiltDownButton = ui.Button("Tilt -", clicked_fn=TiltDown_Click)  

                with ui.HStack():
                    rollRightButton = ui.Button("Roll R", clicked_fn=RollRight_Click)
                    panRightButton = ui.Button("Pan R", clicked_fn=PanRight_Click)
                    tiltUpButton = ui.Button("Tilt +", clicked_fn=TiltUp_Click)                     
                
                # add an IntSlider for translate Strength
                ui.Label("Camera Movement Amount")
                self._movementSlider = ui.IntSlider(min = 10, max = 1000, step=10)                
                self._movementSlider.model.set_value(self._MovementScalar)
                # self._MovementValue = 100
                self._movementSlider.model.add_value_changed_fn(self._on_movement_changed)

                with ui.HStack():
                    upButton = ui.Button("Up", clicked_fn=Up_Click)
                    forwardButton = ui.Button("Forward", clicked_fn=Forward_Click)
                    downButton = ui.Button("Down", clicked_fn=Down_Click)

                with ui.HStack():
                    leftButton = ui.Button("Left", clicked_fn=Left_Click)
                    backButton = ui.Button("Back", clicked_fn=Back_Click)
                    rightButton = ui.Button("Right", clicked_fn=Right_Click)

                with ui.HStack():
                    xAxisButtonUp = ui.Button("X +", clicked_fn=XAxisUp_Click)
                    yAxisButtonUp = ui.Button("Y +", clicked_fn=YAxisUp_Click)
                    zAxisButtonUp = ui.Button("Z +", clicked_fn=ZAxisUp_Click)

                with ui.HStack():
                    xAxisButtonDown = ui.Button("X -", clicked_fn=XAxisDown_Click)
                    yAxisButtonDown = ui.Button("Y -", clicked_fn=YAxisDown_Click)
                    zAxisButtonDown = ui.Button("Z -", clicked_fn=ZAxisDown_Click)

                # with ui.VStack():
                self._label_status_line_1 = ui.Label("")
                self._label_status_line_2 = ui.Label("")
                self._label_buttons = ui.Label("")
                self._label_connected = ui.Label("")
                self._label_debug = ui.Label("")

                # with ui.HStack():
                #     ui.Button("Move", clicked_fn=self.on_click)

        # Note1: It is possible to have multiple 3D mice connected.
        # See: https://github.com/johnhw/pyspacenavigator/blob/master/spacenavigator.py

        self._nav1 = spacenavigator.open(callback=self.on_spacemouse,button_callback=self.on_spacemouse_buttons, DeviceNumber=0)
        self._nav2 = spacenavigator.open(callback=self.on_spacemouse,button_callback=self.on_spacemouse_buttons, DeviceNumber=1)

        if self._nav1 or self._nav2:
            if self._nav1.connected or self._nav2.connected:
                self._label_connected.text = "Connected"
            else:
                self._label_connected.text = "Not Connected"
        else:
            self._label_connected.text = "No spacemouse detected"

    def on_click(self):
        current_time = datetime.now()
        if self.__previous_time:
            if current_time - self.__previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return

        self.__previous_time = current_time
        state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
            **{"t": 255.0, "x": 30.0, "y": 30.0, "z": 30.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
        )
        self.update_state(state)

    def on_spacemouse(self, state: spacenavigator.SpaceNavigator):
        if self.__previous_state == state:
            return
        self.__previous_state = state

        current_time = datetime.now()
        if self.__previous_time:
            if current_time - self.__previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return
            
        self.__previous_time = current_time
        self.update_state(state)

    def on_spacemouse_buttons(self, state: spacenavigator.SpaceNavigator, buttons: spacenavigator.ButtonState):
        current_time = datetime.now()
        if self.__previous_time:
            if current_time - self.__previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return
            
        self.__previous_time = current_time
        self.update_state(state)

    def get_projection_matrix(self, fov, aspect_ratio, z_near, z_far) -> omni.ui.scene.Matrix44:
        """
        Calculate the camera projection matrix.

        Args:
            fov (float): Field of View (in radians)
            aspect_ratio (float): Image aspect ratio (Width / Height)
            z_near (float): distance to near clipping plane
            z_far (float): distance to far clipping plane

        Returns:
            (UsdGeom.Matrix4d): Flattened `(4, 4)` view projection matrix
        """
        a = -1.0 / math.tan(fov / 2)
        b = -a * aspect_ratio
        c = z_far / (z_far - z_near)
        d = z_near * z_far / (z_far - z_near)
        return omni.ui.scene.Matrix44(
            a, 0.0, 0.0, 0.0,
            0.0, b, 0.0, 0.0,
            0.0, 0.0, c, 1.0,
            0.0, 0.0, d, 0.0
        )

    def gfmatrix_to_matrix44(self, matrix: Gf.Matrix4d) -> omni.ui.scene.Matrix44:
        """
        A helper method to convert Gf.Matrix4d to omni.ui.scene.Matrix44

        Args:
            matrix (Gf.Matrix): Input matrix

        Returns:
            UsdGeom.Matrix4d: Output matrix
        """
        # convert the matrix by hand
        # USING LIST COMPREHENSION IS VERY SLOW (e.g. return [item for sublist
        # in matrix for item in sublist]), which takes around 10ms.
        matrix44 = omni.ui.scene.Matrix44(
            matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3],
            matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3],
            matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3],
            matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]
        )
        return matrix44
    
    def gfmatrix_to_array(self, matrix: Gf.Matrix4d) -> list:
        """
        A helper method to convert Gf.Matrix4d to omni.ui.scene.Matrix44

        Args:
            matrix (Gf.Matrix): Input matrix

        Returns:
            UsdGeom.Matrix4d: Output matrix
        """
        # flatten the matrix by hand
        # USING LIST COMPREHENSION IS VERY SLOW (e.g. return [item for sublist
        # in matrix for item in sublist]), which takes around 10ms.
        return (
            matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3],
            matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3],
            matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3],
            matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]
        )

    def decompose_matrix(self, mat: Gf.Matrix4d):
        reversed_ident_matrix = reversed(Gf.Matrix3d())
        translate: Gf.Vec3d = mat.ExtractTranslation()
        scale: Gf.Vec3d = Gf.Vec3d(*(v.GetLength() for v in mat.ExtractRotationMatrix()))
        mat.Orthonormalize()
        rotate: Gf.Vec3d = Gf.Vec3d(*reversed(mat.ExtractRotation().Decompose(*reversed_ident_matrix)))
        return translate, rotate, scale

    def update_translate(self, state, cam_state, update_world_position=False):
        ## On spacemouse, by default x is left(-)/right(+), y is forward(+)/backwards(-), z is up(+)/down(-)
        ## In Omniverse, -z is always camera forwards
        cam_state = VpCamera()
        # Get the current position and target
        cam_pos = cam_state.position_world
        cam_target = cam_state.target_world

        # Create the vector transform - set to state * d
        transform = Gf.Vec3d(
            round(state.x * self._MovementScalar, 1),
            round(state.z * self._MovementScalar, 1),
            round(-state.y * self._MovementScalar, 1)
        )

        if not update_world_position:
            # compute world transform from local
            world_translation = cam_state.usd_camera.ComputeLocalToWorldTransform(Usd.TimeCode.Default()).Transform(transform)
            omni.kit.commands.execute('ChangeProperty',
                                    prop_path=Sdf.Path(self.__camera_path+'.xformOp:translate'),
                                    value=world_translation,
                                    prev=cam_pos)
        else:
            world_translation = transform
            cam_pos = cam_pos + world_translation
            cam_target = cam_target + world_translation

            # Update the world
            cam_state.set_position_world(cam_pos, False)
            cam_state.set_target_world(cam_target, False)
        return transform

            
    def update_rotate(self, state, cam_state, world=False):
            # Get the local transformation - I think we should be using ComputeLocalToWorldTransform rather than GetLocalTransformation & decompose_matrix

            yawsign = 1

            local_transformation: Gf.Matrix4d = cam_state.usd_camera.GetLocalTransformation()
            # translation: Gf.Vec3d = local_transformation.ExtractTranslation()
            # rotation: Gf.Rotation = local_transformation.ExtractRotation()

            decomposed_transform = self.decompose_matrix(local_transformation)

            rotationX = round(decomposed_transform[1][0], 1)
            rotationY = round(decomposed_transform[1][1], 1)
            rotationZ = round(decomposed_transform[1][2], 1)

            # Attempt to hack around issue with going beyond yaw (pan) -90 or +90
            # if( yawsign == 1 and rotationX == -180.0 ):
            #     yawsign = -1
            # elif( yawsign == 1 and rotationX == 180.0 ):
            #     yawsign = -1

            prev_rotation = Gf.Vec3f(rotationX, rotationY, rotationZ)

            new_rotationX = round(rotationX - state.pitch * self._RotationScalar, 1)
            new_rotationY = round(rotationY - state.yaw * self._RotationScalar * yawsign, 1)
            new_rotationZ = round(rotationZ + state.roll * self._RotationScalar, 1)
            alt_local_rotation = Gf.Vec3d(new_rotationX, new_rotationY, new_rotationZ)

            if DEBUG:
                new_rotation = Gf.Rotation(Gf.Vec3d(1, 0, 0), new_rotationX) * \
                                Gf.Rotation(Gf.Vec3d(0, 1, 0), new_rotationY) * \
                                Gf.Rotation(Gf.Vec3d(0, 0, -1), new_rotationZ)
                rotation_transform = Gf.Matrix4d().SetRotate(new_rotation)

                reversed_ident_mtx = reversed(Gf.Matrix3d())
                rotation_transform.Orthonormalize()
                local_rotation = Gf.Vec3d(*reversed(rotation_transform.ExtractRotation().Decompose(*reversed_ident_mtx)))

                self._label_debug.text = f"{new_rotationX:.03f} | {new_rotationY:.03f} | {new_rotationZ:.03f} | {yawsign}"
                self._label_debug.text = self._label_debug.text + '\n' + f"{local_rotation[0]:.03f} | {local_rotation[1]:.03f} | {local_rotation[2]:.03f}"

            world_rotation = alt_local_rotation

            # Update the world
            omni.kit.commands.execute('ChangeProperty',
                                        prop_path=Sdf.Path(self.__camera_path+'.xformOp:rotateXYZ'),
                                        value=world_rotation,
                                        prev=prev_rotation)

    def update_state(self, state: spacenavigator.SpaceNavigator, world=False):
        status_line_1 = f"{state.x:.03f}, {state.y:.03f}, {state.z:.03f}"
        status_line_2 = f"roll: {state.roll:.03f}, tilt: {state.pitch:.03f}, pan: {state.yaw:.03f}"
        # Note1: The number of buttons varies with the type of 3DConnexion product we have
        # Note2: The mappings of buttons is user-configurable so not guaranteed order - we have to account for this
        buttons = f"buttons: {state.buttons}"
        self._label_status_line_1.text = status_line_1
        self._label_status_line_2.text = status_line_2
        self._label_buttons.text = buttons
        self._label_connected.text = f"{state.t}"

        if (
            state.x != 0
            or state.y != 0
            or state.z != 0
            or state.roll != 0
            or state.pitch !=0
            or state.yaw != 0
        ):
            
            ## On spacemouse, by default x is left(-)/right(+), y is forward(+)/backwards(-), z is up(+)/down(-)
            ## In Omniverse, -z is always camera forwards
            cam_state = VpCamera()

            # Update position
            self.update_translate(state, cam_state, world)

            # Now calculate the rotation
            self.update_rotate(state, cam_state, world)

                
    def _on_movement_changed(self, model: ui.SimpleIntModel):
        self._MovementScalar = model.get_value_as_int()
        self._label_debug.text = "Camera movement value = " + str(self._MovementScalar)

    def _onrotation_value_changed(self, model: ui.SimpleIntModel):
        self._RotationValue = model.get_value_as_int()
        self._label_debug.text = "Camera rotation value = " + str(self._RotationValue)

    def on_shutdown(self):
        if self._nav1 is not None:
            self._nav1.close()
            self._nav1.callback = None
            self._nav1.button_callback = None
        if self._nav2 is not None:
            self._nav2.close()
            self._nav2.callback = None
            self._nav2.button_callback = None
        self._nav1 = None
        self._nav2 = None
        self.__previous_time = None
        if self._label_status_line_1 is not None:
            self._label_status_line_1.text = ""
        if self._label_status_line_2 is not None:
            self._label_status_line_2.text = ""
        if self._label_buttons is not None:
            self._label_buttons.text = ""
        if self._label_connected is not None:
            self._label_connected.text = "Not connected"
        self._window = None

        self._active_viewport_window = None
        self._ext_id = None
        print("[coderage.io.spacemouse] coderage io spacemouse shutdown")
