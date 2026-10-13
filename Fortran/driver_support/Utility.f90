FUNCTION find_string_in_array(target_string, string_array) RESULT(index)

  ! Arguments:
  !   target_string (CHARACTER, INTENT(IN)) :: The string to search for.
  !   string_array  (CHARACTER, DIMENSION(:), INTENT(IN)) :: The array of strings to search within.
  !
  ! Returns:
  !   index (INTEGER) :: The 1-based index of the first occurrence of the
  !                      target_string in string_array. Returns 0 if not found.

  CHARACTER(LEN=*), INTENT(IN) :: target_string
  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: string_array
  INTEGER :: index
  INTEGER :: i

  index = 0 ! Initialize index to 0 (not found)

  DO i = 1, SIZE(string_array)
    IF (TRIM(string_array(i)) == TRIM(target_string)) THEN
      index = i
      EXIT ! Exit the loop once the string is found
    END IF
  END DO

END FUNCTION find_string_in_array