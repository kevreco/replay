#include "keys.h"

namespace Re {

Keys::State Keys::states[Keys::KeyArraySize] = {};

int Keys::map[Key_Count] = {};

void Keys::Initialize() {


    for (int i = 0; i < KeyArraySize; ++i) {
        states[i] = Released;
    }

    for (int i = 0; i < Key_Count; ++i) {
        map[i] = -1;
    }
}
void Keys::FrameBegin() {


}

void Keys::FrameEnd() {

    for (int i = 0; i < KeyArraySize; ++i) {

        State& s = states[i];

         if(s == State::OnPress) {
             s = State::Down;
         }
         else if(s == State::OnRelease) {
             s = State::Released;
         }
    }
}


bool Keys::IsDown(Key key) {

    const int mapped_key = map[key];
    const State state = states[mapped_key];

    return (state == OnPress || state == Down);
}

bool Keys::OnPressed(Key key) {

    const int mapped_key = map[key];
    const State state = states[mapped_key];

    return (state == OnPress);
}


} // namespace Re
