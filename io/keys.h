#ifndef RE_KEYS_H
#define RE_KEYS_H

namespace Re {

// The aim of this enum ReKey is to be used like Key::Pressed(KEY_Escape)
// Use the entire set of key could be overkill, we only want the main ones.
// Feel free to remove any key unused in your app.
// Those keys should be registered in Key::Map on the initialization (before the "render loop")

enum Modifier {
    Mod_Ctrl,
    Mod_Alt,
    Mod_Shift
};

enum Key {

    Key_Escape,
    Key_Tab,

    // Arrow
    Key_Left, Key_Up, Key_Right, Key_Down,

    // Numbers
    Key_0, Key_1, Key_2, Key_3,
    Key_4, Key_5, Key_6, Key_7,
    Key_8, Key_9,

    // Letters
    Key_A, Key_B, Key_C, Key_D,
    Key_E, Key_F, Key_G, Key_H,
    Key_I, Key_J, Key_K, Key_L,
    Key_M, Key_N, Key_O, Key_P,
    Key_Q, Key_R, Key_S, Key_T,
    Key_U,Key_V, Key_W, Key_X,
    Key_Y, Key_Z,

    Key_Count

}; // enum Key


// Does not use modifiers, they are not really cross-platform
struct Keys
{

    enum State {

        Released,
        Down,
        OnRelease,
        OnPress
    };

    static const int KeyArraySize = 512;

    static void Initialize();

    static void FrameBegin();
    static void FrameEnd();

    static bool IsDown(Key key);

    // True only for the first frame where the key has been pressed
    static bool OnPressed(Key key);

    // OnPress, OnRelease etc
    static State states[KeyArraySize];

    // Map from Re::Key enum to system numeric value;
    // Should be initialized at application start.
    static int map[Key_Count];

}; // struct Keys

} // namespace Re

/* Example, using GLFW

-  First a function must register all the desired keys at the application initialization

void InitializeKeyboard() {

    Re::Keys::map[Re::Key_Tab] = GLFW_KEY_TAB;
    Re::Keys::map[Re::Key_Escape] = GLFW_KEY_ESCAPE;

    Re::Keys::map[Re::Key_Z] = GLFW_KEY_Z;
    Re::Keys::map[Re::Key_Q] = GLFW_KEY_Q;
    Re::Keys::map[Re::Key_S] = GLFW_KEY_S;
    Re::Keys::map[Re::Key_D] = GLFW_KEY_D;

    Re::Keys::map[Re::Key_Up] = GLFW_KEY_UP;
    Re::Keys::map[Re::Key_Right] = GLFW_KEY_RIGHT;
    Re::Keys::map[Re::Key_Down] = GLFW_KEY_DOWN;
    Re::Keys::map[Re::Key_Left] = GLFW_KEY_LEFT;
}


- Then, a function must setup the pressed/released key on the begining of the frame

void KeyCallback(GLFWwindow* w, int key, int scancode, int action, int mods) {

    // Setup Re Keys

    if (action == GLFW_PRESS) {
        Keys::State& s = Keys::states[key];

        if (s != Keys::Down) {
            s = Keys::OnPress;
        }

    }

    else if (action == GLFW_RELEASE) {
        Keys::State& s = Keys::states[key];
        if (s != Keys::Released) {
            s = Keys::OnRelease;
        }
    }
}

- Finally Keys::FrameEnd() must be called at the frame end

Keys::FrameEnd();

*/

#endif // RE_KEYS_H
