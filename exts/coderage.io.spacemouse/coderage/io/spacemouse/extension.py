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
from omni.kit.widget.viewport.api import ViewportAPI
from omni.kit.viewport.window import ViewportWindow
from pxr import Usd, UsdGeom, Gf
if platform.system() == 'Windows':
    omni.kit.pipapi.install("pywinusb")
# import pywinusb
import spacenavigator

UPDATE_TIME_MILLIS = 500


class CoderageIoSpacemouseExtension(omni.ext.IExt):
    def __init__(self):
        self._count = 0
        self.previous_time = None
        self.previous_state = None

    def on_startup(self, ext_id):
        print("[coderage.io.spacemouse] coderage io spacemouse startup")

        self._window = ui.Window("Spacemouse debug", width=300, height=300)
        with self._window.frame:
            with ui.VStack():
                self._label_msg = ui.Label("")
                self._label_msg2 = ui.Label("")
                self._label_buttons = ui.Label("")
                self._label_connected = ui.Label("")

                with ui.HStack():
                    ui.Button("Move", clicked_fn=self.on_click)

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
        if self.previous_time:
            if current_time - self.previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return

        self.previous_time = current_time
        state: spacenavigator.SpaceNavigator = spacenavigator.SpaceNavigator(
            **{"t": 255.0, "x": 30.0, "y": 30.0, "z": 30.0, "roll": 0.0, "pitch": 0.0, "yaw": 0.0, "buttons": [0,0]}
        )
        self.update_state(state)

    def on_spacemouse(self, state: spacenavigator.SpaceNavigator):
        if self.previous_state == state:
            return
        self.previous_state = state

        current_time = datetime.now()
        if self.previous_time:
            if current_time - self.previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return
            
        self.previous_time = current_time
        self.update_state(state)

    def on_spacemouse_buttons(self, state: spacenavigator.SpaceNavigator, buttons: spacenavigator.ButtonState):
        current_time = datetime.now()
        if self.previous_time:
            if current_time - self.previous_time < timedelta(milliseconds=UPDATE_TIME_MILLIS):
                return
            
        self.previous_time = current_time
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

    def update_state(self, state: spacenavigator.SpaceNavigator):
        msg = f"{state.x:.03f}, {state.y:.03f}, {state.z:.03f}"
        msg2 = f"roll: {state.roll:.03f}, pitch: {state.pitch:.03f}, yaw: {state.yaw:.03f}"
        # Note1: The number of buttons varies with the type of 3DConnexion product we have
        # Note2: The mappings of buttons is user-configurable so not guaranteed order - we have to account for this
        buttons = f"buttons: {state.buttons}"
        self._label_msg.text = msg
        self._label_msg2.text = msg2
        self._label_buttons.text = buttons
        self._label_connected.text = f"{state.t}"

        ctx = omni.usd.get_context()
        stage: Usd.Stage = ctx.get_stage()

        active_viewport_window: ViewportWindow = ViewportWindow.active_window
        active_viewport_api: ViewportAPI = active_viewport_window.viewport_api

        # active_camera_path = viewport_api.get_active_camera()
        if active_viewport_window:
            viewport_time = active_viewport_api.time
            active_camera_path = active_viewport_api.camera_path

        active_camera_prim: Usd.Prim = stage.GetPrimAtPath(active_camera_path)

        # usd_camera = UsdGeom.Camera(active_camera_prim)
        # self._gf_camera: Gf.Camera = usd_camera.GetCamera()
        
        usd_camera = UsdGeom.Xformable(active_camera_prim)

        projection = self.gfmatrix_to_array(usd_camera.GetLocalTransformation())

        if (
            state.x != 0
            or state.y != 0
            or state.z != 0
        ):
            # Based on code from
            # https://github.com/mati-nvidia/developer-office-hours/blob/main/exts/maticodes.doh_2023_04_14/scripts/move_prim_forward.py
            local_transformation: Gf.Matrix4d = usd_camera.GetLocalTransformation()
            # Apply the local matrix to the start and end points of the camera's default forward vector (-Z)
            a: Gf.Vec4d = Gf.Vec4d(0,0,0,1) * local_transformation
            b: Gf.Vec4d = Gf.Vec4d(0,0,-10,1) * local_transformation
            # Get the vector between those two points to get the camera's current forward vector
            cam_fwd_vec = b-a
            # Convert to Vec3 and then normalize to get unit vector
            cam_fwd_unit_vec = Gf.Vec3d(cam_fwd_vec[:3]).GetNormalized()
            # Multiply the forward direction vector with how far forward you want to move
            forward_step = cam_fwd_unit_vec * 100
            # Create a new matrix with the translation that you want to perform
            offset_mat = Gf.Matrix4d()
            offset_mat.SetTranslate(forward_step)
            # Apply the translation to the current local transform
            new_transform: Gf.Matrix4d = local_transformation * offset_mat
            # Extract the new translation
            translate: Gf.Vec3d = new_transform.ExtractTranslation()

            # Update the attribute - this next line hangs
            active_camera_prim.GetAttribute("xformOp:translate").Set(translate)


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
        self.previous_time = None
        if self._label_msg is not None:
            self._label_msg.text = ""
        if self._label_msg2 is not None:
            self._label_msg2.text = ""
        if self._label_buttons is not None:
            self._label_buttons.text = ""
        if self._label_connected is not None:
            self._label_connected.text = "Not connected"
        self._window = None

        self._active_viewport_window = None
        self._ext_id = None
        print("[coderage.io.spacemouse] coderage io spacemouse shutdown")
