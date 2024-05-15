// Rasperry Pi Pico daemon reading data from GPIO pin(s) and sending it via the serial port

// initialize variables and constants
const int buttonPin1 = 27;    // the number of the pushbutton pin
const int buttonPin2 = 28;    // the number of the pushbutton pin
int buttonState1;             // the current reading from the input pin
int lastButtonState1 = LOW;   // the previous reading from the input pin
int buttonState2;             // the current reading from the input pin
int lastButtonState2 = LOW;   // the previous reading from the input pin

// the following variables are unsigned longs because the time, measured in
// milliseconds, will quickly become a bigger number than can be stored in an int
unsigned long lastDebounceTime1 = 0;  // the last time the output pin was toggled
unsigned long lastDebounceTime2 = 0;  // the last time the output pin was toggled
unsigned long debounceDelay = 50;     // the debounce time; increase if the output flickers

void setup() {
  pinMode(buttonPin1, INPUT);
  pinMode(buttonPin2, INPUT);
  Serial.begin(115200);
}

void loop() {
  // read the state of the switch into a local variable:
  int reading1 = digitalRead(buttonPin1);
  int reading2 = digitalRead(buttonPin2);

  // check to see if you just pressed the button
  // (i.e. the input went from LOW to HIGH), and you've waited long enough
  // since the last press to ignore any noise:

  // if the switch changed, due to noise or pressing:
  if (reading1 != lastButtonState1) {
    // reset the debouncing timer
    lastDebounceTime1 = millis();
  }
  if (reading2 != lastButtonState2) {
    // reset the debouncing timer
    lastDebounceTime2 = millis();
  }

  if ((millis() - lastDebounceTime1) > debounceDelay) {
    // whatever the reading is at, it's been there for longer than the debounce
    // delay, so take it as the actual current state:

    // if the button state has changed:
    if (reading1 != buttonState1) {
      buttonState1 = reading1;

      // send the button state to Serial port
      if (buttonState1 == HIGH) {
        Serial.print(buttonPin1);
        Serial.print(":");
        Serial.println("1");
      } else {
        Serial.print(buttonPin1);
        Serial.print(":");
        Serial.println("0");
      }
    }
  }

  if ((millis() - lastDebounceTime2) > debounceDelay) {
    // whatever the reading is at, it's been there for longer than the debounce
    // delay, so take it as the actual current state:

    // if the button state has changed:
    if (reading2 != buttonState2) {
      buttonState2 = reading2;

      // send the button state to Serial port
      if (buttonState2 == HIGH) {
        Serial.print(buttonPin2);
        Serial.print(":");
        Serial.println("1");
      } else {
        Serial.print(buttonPin2);
        Serial.print(":");
        Serial.println("0");
      }
    }
  }

  // save the reading
  lastButtonState1 = reading1;
  lastButtonState2 = reading2;
}
